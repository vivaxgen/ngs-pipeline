Presenting Sample Read Files
============================

To be able to process the FASTQ read files, NGS-Pipeline needs to access those
read files.
There are several mechanisms on how FASTQ read files can be presented to the
pipeline.
Which mechanism is better is left to the discretion of users, since this will
depend on the circumstances and scenario.

If the FASTQ read files are directly accesible from the working directory
(either locally or remotely through networked or clustered/parallel file
system), symbolic links can be generated inside user's working/analysis
directory, as follow:

* Directly link the source directory containing the read files

  This is the simplest step, using the following command::

    ln -s /PATH/TO/SOURCE/READ_DIR reads

* Generating symbolic links to individual read files

  In this scheme, rather than generating a link to the source directory
  containing read files, individual symbolic link to every read files are
  generated, eg::

    mkdir reads
    ln -sr /PATH/TO/SOURCE/READ_DIR/sample-1_R1.fastq.gz reads/
    ln -sr /PATH/TO/SOURCE/READ_DIR/sample-2_R2.fastq.gz reads/
    ...

  BASH for-loop can be utilized to perform the above steps in a more automatic
  way::

    mkdir reads
    for fn in /PATH/TO/SOURCE/READ_DIR/sample*.fastq.gz; do ln -sr {fn} reads/; done

  Other alternative is to use ``ngs-pl generate-links`` command as follows::

    mkdir reads
    ngs-pl generate-links -o reads /PATH/TO/SOURCE/READ_DIR/sample*.fastq.gz

  The ``generate-links`` command can accept multiple sources from different
  directories as well::

    ngs-pl generate-links -o reads /DIR_1/*.fastq.gz /DIR_2/*.fastq.gz

  Some commands in NGS-Pipeline require that the sample codes are used as the
  filenames of the FASTQ files, eg. ``sample-123_R1.fastq.gz`` (for paired-end
  input) or ``sample-ABC.fastq.gz`` (for single read input).
  A feature from ``generate-links`` is the ``--underscore`` option that
  can be used in cases where the filenames of the read files contain more than
  the actual sample code.
  Some sequencing providers add the well number, or the lane number, or even
  the indexing sequences to be used.
  For example, the provider might named the read files as below::

    ID004289_RVD0_L9V2F_CCAACCGTTC-GCAAGAGATG_L001_R1.fastq.gz

  The filename above contains a run serial number ``L9V2F``, the indexing
  barcode ``CCAACCGTTC-GCAAGAGATG`` and the lane number ``L001`` with the
  actual sample name of just ``ID004289_RVD0``.
  The underscore that separate the sample name and the rest of the information
  is the 4th underscore, counting from the end of filename.
  In this case, by using the following example commands::

    ngs-pl generate-links -o reads --underscore 4 SOURCE_DIR/*.fastq.gz

  the layout will be::

    SOURCE_DIR/ID004289_RVD0_L9V2F_CCAACCGTTC-GCAAGAGATG_L001_R1.fastq.gz

    reads/ID004289_RVD0_R1.fastq.gz -> SOURCE_DIR/ID004289_RVD0_L9V2F_CCAACCGTTC-GCAAGAGATG_L001_R1.fastq.gz

  Counting underscore from the end of the filename will be consistent since
  sample name might contain underscore as well.

If the read files are not directly accessible by NGS-Pipeline, the read files
need to be copied to a directory accessible by NGS-Pipeline.
Although the read files can be copied to the current working directory,
one of the common practice is to copy those read files to a directory that
will hold all FASTQ read files, and then create the necessary symbolic links.


