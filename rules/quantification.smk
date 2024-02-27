## Performs transcript quantification on RNA-seq data
rule salmon:
    input:
        fa=rules.retrieve_transcripts_fa.output.fa,
        bam=rules.run_star.output.bam_transcriptome
    output:
        sf=paths.salmon.sf
    benchmark:
        'benchmark/{sample}_salmon.tab'
    log:
        'log/{sample}_salmon.log'
    conda:
        SOURCEDIR+"/../envs/salmon.yaml"
    params:
        path=PREDIR+"/salmon/{sample}"
    threads: max(1,min(8,NCORES))
    shell:
        '''
          echo "salmon quant -t {input.fa} -l A -a {input.bam} -p {threads} -o {params.path} && mv {params.path}/quant.sf {output.sf}" | tee {log}
          salmon quant -t {input.fa} -l A -a {input.bam} -p {threads} -o {params.path} && mv {params.path}/quant.sf {output.sf} 2>> {log}
         
          ## Export rule env details
          conda env export --no-builds > info/salmon.info 
        '''
