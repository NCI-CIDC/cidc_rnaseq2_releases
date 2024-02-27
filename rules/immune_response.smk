## Run MSIsensor2 to assess microsatellite instability (MSI) in tumor samples
rule msisensor2:
    input:
        bam=rules.run_star.output.bam,
        bai=rules.index_bam.output,
        models=rules.retrieve_immune_refs.output.models,
        tch=rules.retrieve_immune_refs.output.tch
    output:
        output=paths.msisensor2.output,
        dis=paths.msisensor2.dis,
        somatic=paths.msisensor2.somatic,
        txt=paths.msisensor2.txt
    benchmark:
        'benchmark/{sample}_msisensor2.tab'
    log:
        'log/{sample}_msisensor2.log'
    conda:
        SOURCEDIR+"/../envs/msisensor2.yaml"
    threads: max(1,min(8,NCORES))
    params:
        prefix='msisensor2/{sample}_msisensor2'
    shell:
        '''
          echo "msisensor2 msi -M {input.models} -t {input.bam} -o {params.prefix} -b {threads} && cp {output.output} {output.txt}" | tee {log}
          msisensor2 msi -M {input.models} -t {input.bam} -o {params.prefix} -b {threads} && cp {output.output} {output.txt} 2>> {log}
          
          ## Export rule env details
          conda env export --no-builds > info/msisensor2.info
        '''
