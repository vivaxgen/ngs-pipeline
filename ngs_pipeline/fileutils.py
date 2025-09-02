from enum import Enum
import pathlib
from ngs_pipeline import cerr

ReadMode = Enum("ReadMode", ["SINGLETON", "PAIRED_END"])


class ReadFileDict(object):
    """a dictionary with sample as keys, and values of [file1, ...] for singleton
    or [(file_R1, file_R2), ...] for paired-end
    """

    ReadMode = ReadMode

    def __init__(
        self,
        infiles: list,
        underscore: int,
        underscore_prefix: int = 0,
        remove_prefix: str | None = None,
        mode: str | None = None,
        skip_list: list = [],
    ):
        super().__init__()
        self._d = {}
        self.err_files = []
        self.mode = mode
        self.underscore = underscore
        self.underscore_prefix = underscore_prefix
        self.remove_prefix = remove_prefix
        self.skip_list = skip_list
        self.populate_read_files(infiles)

    def keys(self):
        return self._d.keys()

    def __getitem__(self, key):
        return self._d[key]

    def __setitem__(self, key, value):
        self._d[key] = value

    def __contains__(self, key):
        return key in self._d

    def __iter__(self):
        raise NotImplementedError()

    def samples(self):
        return list(sorted(self.keys()))

    def get_read_file(self, wildcards):
        """this function is suitable to be used in input: directive in snake rules
        wildcards must contain "sample" field
        """
        idx = int(wildcards.idx)
        return self[wildcards.sample][idx]

    def get_read_file_as_dict(self, wildcards):
        """this function is suitable to be used in input: directive in snake rules
        wildcards must contain "sample" field
        return a dict of {'read1': , 'read2': }
        usage:
            input:
                unpack(obj.get_read_file_as_dict)
        """
        idx = int(wildcards.idx)
        read_files = self[wildcards.sample][idx]
        if len(read_files) == 1:
            return dict(read0=read_files[0])
        elif len(read_files) > 1:
            d = {}
            for i in range(len(read_files)):
                d["read" + str(i + 1)] = read_files[i]
        else:
            raise ValueError("the sample does not have read files")

        print(d)
        return d

    def get_indexes_with_wildcards(self, wildcards):
        """this function is suitable to be used in input: directive in snake rules
        wildcards must contain "sample" field
        return a list of read indexes for each sample
        """
        return self.get_indexes(wildcards.sample)

    def get_indexes(self, sample):
        """return list of indexes for each sample"""
        return list(range(len(self[sample])))

    def populate_read_files(self, infiles):

        if len(infiles) == 0:
            raise ValueError("No input files given")

        # check the existence of the 1st file to determine that the rest of the files
        # are also exist from shell globbing
        if not pathlib.Path(infiles[0]).exists():
            raise ValueError(
                f"No input files given as the first file is [{infiles[0]}]"
            )

        infiles = sorted(infiles)
        self.mode = check_read_mode(infiles) if not self.mode else self.mode

        err_files = self.err_files
        if self.mode == ReadMode.SINGLETON:

            for infile in infiles:
                if not infile.endswith(".fastq.gz"):
                    err_files.append(
                        f"Input file: {infile} is not a compressed FASTQ file (fastq.gz)."
                    )
                    continue

                path = pathlib.Path(infile)
                if not path.exists():
                    err_files.append(
                        f"Input file: {infile} does not exist. Please check your path."
                    )
                    continue

                sample = get_sample_name(
                    infile, self.underscore, self.underscore_prefix, self.remove_prefix
                )
                if sample in self.skip_list:
                    continue
                if sample not in self:
                    self[sample] = [(infile,)]
                else:
                    self[sample].append((infile,))

        elif self.mode == ReadMode.PAIRED_END:

            if (len(infiles) % 2) != 0:
                raise ValueError(
                    f"Error: reads are paired-end but number of infiles is odd ({len(infiles)})"
                )

            # with paired-end reads, underscore must at least be 1, since filenames have _R1.fastq.gz
            # and _R2.fastq.gz (or _!.fastq.gz and _2.fastq.gz)
            if self.underscore == 0:
                self.underscore = 1

            infiles_1, infiles_2 = infiles[::2], infiles[1::2]
            for infile_1, infile_2 in zip(infiles_1, infiles_2):
                prefix_1 = get_sample_name(
                    infile_1,
                    self.underscore,
                    self.underscore_prefix,
                    self.remove_prefix,
                )
                prefix_2 = get_sample_name(
                    infile_2,
                    self.underscore,
                    self.underscore_prefix,
                    self.remove_prefix,
                )
                if prefix_1 != prefix_2:
                    err_files.append(
                        f"ERROR: unmatch pair [{prefix_1}] <> [{prefix_2}]"
                    )

                if prefix_1 in self.skip_list:
                    continue
                if prefix_1 not in self:
                    self[prefix_1] = [(infile_1, infile_2)]
                else:
                    self[prefix_1].append((infile_1, infile_2))
        else:
            raise ValueError(f"Unknown mode: {self.mode}")


def check_read_mode(infiles: list):
    """check if infiles are paired-rend or singleton reads"""

    if (len(infiles) % 2) != 0:
        cerr(
            f"Number of infiles is not even [{len(infiles)}], assuming singleton reads."
        )
        return ReadMode.SINGLETON

    counters = {
        "1.fastq.gz": 0,
        "2.fastq.gz": 0,
    }
    for infile in infiles:

        if not "_" in infile:
            cerr(
                "Filename does not contain underscore character, hence assumming singleton reads."
            )
            return ReadMode.SINGLETON

        if (suffix := infile.rsplit("_", 1)[-1].removeprefix("R")) in counters:
            counters[suffix] += 1

    if sum(counters.values()) != len(infiles):
        cerr(
            "Number of infiles ending with 1.fastq.gz or 2.fastq.gz is not equal to "
            "total files, assuming singleton reads."
        )
        return ReadMode.SINGLETON

    if counters["1.fastq.gz"] != counters["2.fastq.gz"]:
        cerr(
            "Number of infiles wnding with 1.fastq.gz is not equal with infiles "
            "ending with 2.fastq.gz, assuming singleton reads"
        )
        return ReadMode.SINGLETON

    return ReadMode.PAIRED_END


def get_sample_name(
    filename, underline=0, remove_underline_prefix=0, remove_prefix=None
):
    """
    get sample name after splitting a number of underline character
    and removing prefix
    """
    filename = pathlib.Path(filename)
    # if contain tilde, just split on the tilde
    if "~" in filename.name:
        return filename.name.split("~")[0]
    if remove_prefix:
        filename = pathlib.Path(filename.name.removeprefix(remove_prefix))
    if remove_underline_prefix:
        filename = pathlib.Path(filename.name.split("_", remove_underline_prefix)[-1])
    if underline != 0:
        return filename.name.rsplit("_", underline)[0]
    return filename.name.removesuffix(".fastq.gz")


def create_relative_symlink(dest: pathlib.Path, source: pathlib.Path, force=False):

    if force:
        dest.unlink(missing_ok=True)

    # find common parent directory
    dest_abs = dest.resolve()
    source_abs = source.absolute()

    # create a relative path from the destination file to the source file
    relative_src = source_abs.relative_to(dest_abs.parent, walk_up=True)

    # create a symbolic link from dest_file to the relative source path
    dest_abs.symlink_to(relative_src)

    return

    common_dir = None
    for d1, d2 in zip(reversed(dest_abs.parents), reversed(source_abs.parents)):
        if d1 == d2:
            common_dir = d1
        else:
            break

    # find relative path from common dir
    dest_from_common = dest_abs.relative_to(common_dir)
    source_from_common = source_abs.relative_to(common_dir)

    # generate symbolic link
    rel_link = pathlib.Path("/".join([".."] * (len(dest_from_common.parents) - 1)))
    rel_link = rel_link / source_from_common

    dest_abs.symlink_to(rel_link)


def _make_symlink_for_sample(
    sample: str,
    path: str,
    outdir: pathlib.Path,
    ext: str,
    *,
    use_absolute: bool = False,
    overwrite_if_exists: bool = False,
):

    # check that path exists
    path = pathlib.Path(path)
    if not path.exists():
        raise ValueError(f"File {path} does not exist!")

    target_path = outdir / (sample + ext)
    if target_path.exists():
        if not overwrite_if_exists:
            raise ValueError(f"File {target_path} already exists!")

        # check if the existing link is identical with the link
        if path.absolute().resolve() == target_path.absolute().resolve():
            # link and path is identical, just skip
            return
        target_path.unlink()

    # find commmon denominator directory of path and target path
    if use_absolute:
        target_path.symlink_to(path.resolve())
    else:
        create_relative_symlink(target_path, path)


def make_sample_symlink(
    sample: str,
    file_item: tuple,
    outdir: str | pathlib.Path,
    mode: ReadMode = ReadMode.PAIRED_END,
    index: int = -1,
    *,
    use_absolute: bool = False,
):

    if not type(outdir) == pathlib.Path:
        outdir = pathlib.Path(outdir)

    if not outdir.exists():
        raise ValueError(f"Output directory {outdir} does not exist!")

    if type(file_item) != tuple:
        raise ValueError(f"File {file_item} for sample {sample} is not a tuple")

    sample_code = f"{sample}~{index}" if index >= 0 else sample

    if mode == ReadMode.PAIRED_END:
        if len(file_item) != 2:
            raise ValueError(f"File {file_item} for sample {sample} is paired-end")
        for pair_idx, path in enumerate(file_item, 1):
            _make_symlink_for_sample(
                sample_code,
                path,
                outdir,
                f"_R{pair_idx}.fastq.gz",
                use_absolute=use_absolute,
            )

    else:
        _make_symlink_for_sample(
            sample_code, file_item[0], outdir, f".fastq.gz", use_absolute=use_absolute
        )


# EOF
