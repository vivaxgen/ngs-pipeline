

rule gen_report:
    threads: 1
    input:
        vcf = '{pfx}/vcfs/variants.vcf.gz'
    output:
        tsv = '{pfx}/genetic_report.tsv'
    params:
        extra_flags = config.get('generate_variant_report_extra_flags', '')
    shell:
        "ngs-pl generate-variant-report {params.extra_flags} -o {output.tsv} "
        "--infofile {variant_info} {input.vcf}"


rule merge_report:
    localrule: True
    input:
        expand(f'{outdir}/{{sample}}/genetic_report.tsv', sample=read_files.samples()))
    output:
        f'{outdir}/merged_report.tsv'
    run:
        import pandas as pd

        dfs = []
        for infile in input:
            try:
                df = pd.read_table(infile)
                dfs.append(df)
            except:
                raise RuntimeError(f'ERROR parsing {infile}')

        concatenated_df = pd.concat(dfs)
        concatenated_df.to_csv(output[0], index=False, sep='\t')


rule merge_report_xlsx:
    localrule: True
    input:
        f'{outdir}/merged_report.tsv'
    output:
        f'{outdir}/merged_report.xlsx'
    run:
        import pandas as pd

        df = pd.read_table(input[0])
        df.to_excel(output[0], index=False, sheet_name="Merged_Report")


# EOF
