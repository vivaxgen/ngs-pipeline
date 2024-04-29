# msf_panel_genreport.smk - ngs-pipeline rules
# [https://github.com/vivaxgen/ngs-pipeline]

__copyright__ = "(C) 2023, Hidayat Trimarsanto <trimarsanto@gmail.com>"
__license__ = "MIT"


rule gen_report:
    threads: 1
    input:
        vcf = '{pfx}/vcfs/variants.vcf.gz'
    output:
        tsv = '{pfx}/genetic_report.tsv'
    params:
        extra_flags = config.get('generate_variant_report_extra_flags', '')
        min_var_qual = config.get('min_variant_qual', 30)
        min_depth = config.get('min_depth', 10)
    shell:
        "ngs-pl generate-variant-report --mindepth {params.min_depth} --min-var-qual {params.min_var_qual} {params.extra_flags} -o {output.tsv} "
        "--infofile {variant_info} {input.vcf}"


rule merge_report:
    localrule: True
    input:
        expand(f'{outdir}/{{sample}}/genetic_report.tsv', sample=read_files.samples())
    output:
        f'{outdir}/merged_genetic_report.tsv'
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
        f'{outdir}/merged_genetic_report.tsv'
    output:
        f'{outdir}/merged_genetic_report.xlsx'
    run:
        import pandas as pd

        df = pd.read_table(input[0])
        df.to_excel(output[0], index=False, sheet_name="Merged_Report")


# EOF
