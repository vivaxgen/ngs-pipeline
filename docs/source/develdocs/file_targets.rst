File Targets
============

The `<sp>` prefix is the Snakemake pathvars for <sample prefix>.
It `<sp> = ""` for single-sample flow, and `<sp> = "{outdir}/samples/{sample}/"` for multi-sample flow.

.. list-table:: Target (Input/Output) files for the rules
   :widths: 50 50
   :header-rows: 1

   * - No.
     - Target Files
     - Remarks
   * - 1.
     - ::
          <sp>reads/raw-{idx}_R1.fastq.gz
          <sp>reads/raw-{idx}_R2.fastq.gz

     - The symbolic link to the raw reads files, which is the input for the trimming rules.

   * - 2.
     - ::
          <sp>trimmed-reads/trimmed-{idx}_R1.fastq.gz
          <sp>trimmed-reads/trimmed-{idx}_R2.fastq.gz

     - The output of the trimming rules, which is the input for the mapping rules.

   * - 3.
     - ::

          <sp>maps/mapped-{idx}.bam

     - The output of the mapping rules, which is the input for the variant calling rules.
       This file contains either all the reads except those whose at least one of the pair is mapped to contaminant sequences,
       or all the reads whose at least one of the pair is mapped to the target reference genome.
       The BAM file is name-sorted and fix-mated.
       This is a temporary file, to keep it set `keep_paired_bam = True` in the configuration file.

   * - 4.
     - ::

          <sp>maps/mapped-filtered-{idx}.bam
     - This is the BAM file that has been filtered based on criteria set in the configuration file.
       Filtering that can be performed includes: removing unmapped reads, remove secondary/supplementary/trans maps, and orientation based filtering (FF, FR, RF, RR).
       See the command `filter-reads-orientation` for more details.
       This BAM file is coordinate-sorted.

   * - 5.
     - ::

          <sp>maps/mapped-dedup-{idx}.bam

     - The deduplicated BAM file, which is suitable for non-amplicon (eg. WGS) data set.

   * - 6.
     - ::

          <sp>maps/mapped-final-{idx}.bam

     - This is the final BAM file before being merged (from multiple runs of the same sample).
       This would be a hardlink to either filtered [4] BAM file (`deduplicate = False`) or deduplicated (5) BAM file (`deduplicate = True`).
    
   * - 7.
     - ::

          <sp>maps/mapped-final.bam

     - The merged BAM file, which is the input for the variant calling rules.
       This is a temporary file, to keep it set `keep_final_bam = True` in the configuration file.

    


