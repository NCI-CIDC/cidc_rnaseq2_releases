## Run msisensor
rule msisensor2_tumor_call:
    input:
        bam=rules.run_bwa.output,
        bai=rules.sample_bam.output.index
    output:
        msisensor="msisensor2/{sample}_msisensor2",
        msisensor_txt="msisensor2/{sample}_msisensor2.txt",
        msisensor_dis="msisensor2/{sample}_msisensor2_dis",
        msisensor_somatic="msisensor2/{sample}_msisensor2_somatic"
    benchmark:
        'benchmark/{sample}_msisensor.tab'
    log:
        'log/{sample}_msisensor.log'
    conda:
        SOURCEDIR+"/../envs/immune_response.yaml"
    params:
        sample='{sample}'
        prefix='msisensor2/{sample}_msisensor2',
        msisensor2_models= paths.immune_response.msisensor2_models
    shell:
        '''
          msisensor2 msi -M {params.msisensor2_models} -t {input.bam} -o {params.prefix}
          cp {output.msisensor} {output.msisensor_txt}
        '''
