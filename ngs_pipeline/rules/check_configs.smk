
outfile = config['outfile']

rule all:
    input:
        f'{outfile}'

rule write_configs:
    localrule: True
    output:
        f'{outfile}'
    params:
        config = config
    run:

        import yaml

        del params.config['outfile']
        with open(output[0], 'w') as f_out:
            yaml.dump(params.config, f_out)

# EOF
