
N = int(config.get('N', 4))

rule all:
    input:
        *[ f'test-cluster-{idx}' for idx in range(N)]


rule send_to_cluster:
    localrule: False
    output:
        'test-cluster-{idx}'
    shell:
        'hostname > {output}'


# EOF
