
from ngs_pipeline import cerr

rule map_stats:
    threads: 1
    input:
        "{pfx}/maps/{filename}.bam"
    output:
        "{pfx}/logs/{filename}.stats.txt"
    shell:
        "samtools stats {input} > {output}"


rule map_depth:
    # this rule creates stats and depths of any .bam file using samtools stats
    threads: 1
    input:
        '{pfx}/maps/{filename}.bam'
    output:
        '{pfx}/logs/{filename}.depths.txt.gz'
    shell:
        'samtools depth {input} | gzip > {output}'


rule collect_stats:
    threads: 1
    input:
        trims = expand('{{pfx}}/{{sample}}/logs/trimming_stat-{idx}.json', idx=read_files.get_indexes_with_wildcards),
        maps = expand('{{pfx}}/{{sample}}/logs/{{sample}}-{idx}.stats.txt', idx=read_files.get_indexes_with_wildcards),
        filtered = expand('{{pfx}}/{{sample}}/logs/mapped-filtered-{idx}.stats.txt', idx=read_files.get_indexes_with_wildcards),
        dedups = expand('{{pfx}}/{{sample}}/logs/mapped-dedup-{idx}.stats.txt', idx=read_files.get_indexes_with_wildcards) if deduplicate else [],
        finals = expand('{{pfx}}/{{sample}}/logs/mapped-final-{idx}.stats.txt', idx=read_files.get_indexes_with_wildcards),
        depths = expand('{{pfx}}/{{sample}}/logs/mapped-final-{idx}.depths.txt.gz', idx=read_files.get_indexes_with_wildcards),
    params:
        trimmed = lambda wildcards, input: '--trimmed ' + ' --trimmed '.join(input.trims),
        mapped = lambda wildcards, input: '--mapped ' + ' --mapped '.join(input.maps),
        deduped = (lambda wildcards, input: '--dedup ' + ' --dedup '.join(input.dedups)) if deduplicate else '',
        finaled = lambda wildcards, input: '--final ' + ' --final '.join(input.finals),
        depthed = lambda wildcards, input: '--depth ' + ' --depth '.join(input.depths),
    output:
        stat = '{pfx}/{sample}/logs/stats.tsv'
    shell:
        'ngs-pl calculate-stats -o {output} --mindepth {min_depth} '
        '{params.trimmed} {params.mapped} {params.deduped} {params.finaled} {params.depthed} {wildcards.sample}'


rule gather_stats_sample:
    threads: 1
    input:
        expand(f'{outdir}/samples/{{sample}}/logs/stats.tsv', sample=read_files.samples()),
    output:
        stat = f'{outdir}/stats.tsv'
    run:

        # import heavy modules here if required
        from pathlib import Path
        import pandas as pd

        cerr('Gathering stats from pipeline results')

        dfs = []

        for statfile in input:
                dfs.append(pd.read_table(statfile, sep='\t'))

        # concat all
        cerr(f'[Gathering stats file from {len(dfs)} directories]')
        all_df = pd.concat(dfs)

        all_df.to_csv(output.stat, index=False, sep='\t')


# EOF
