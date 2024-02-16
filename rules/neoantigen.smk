## Extracts reads mapped to chromosome 6 and any HLA decoys or chromosome 6 alternates with arcasHLA
rule arcashla_extract:
    input:
        tch=rules.retrieve_imgthla_db.output.tch,
        bam=rules.run_star.output.bam
    output:
        fq1=paths.arcashla.fq1,
        fq2=paths.arcashla.fq2
    benchmark:
        'benchmark/{sample}_arcashla_extract.tab'
    log:
        'log/{sample}_arcashla_extract.log'
    conda:
        SOURCEDIR+"/../envs/arcashla.yaml"
    params:
        path=PREDIR+"/arcashla/{sample}"
    threads: max(1,min(8,NCORES))
    shell:
        '''
          echo "arcasHLA extract {input.bam} -t {threads} -v -o {params.path}" | tee {log}
          arcasHLA extract {input.bam} -t {threads} -v -o {params.path} 2>> {log}
        '''

## Genotypes HLA alleles from extracted reads (no partial alleles) with arcasHLA
rule aracshla_genotype:
    input:
        fq1=rules.arcashla_extract.output.fq1,
        fq2=rules.arcashla_extract.output.fq2
    output:
        p=paths.arcashla.p,
        json_genes=paths.arcashla.json_genes,
        json_genotype=paths.arcashla.json_genotype,
        log=paths.arcashla.log
    benchmark:
        'benchmark/{sample}_arcashla_genotype.tab'
    log:
        'log/{sample}_arcashla_genotype.log'
    conda:
        SOURCEDIR+"/../envs/arcashla.yaml"
    params:
        path=PREDIR+"/arcashla/{sample}"
    threads: max(1,min(8,NCORES))
    shell:
        '''
          echo "arcasHLA genotype {input.fq1} {input.fq2} -g A,B,C,DQA1,DQB1,DRB1 -t {threads} -v -o {params.path}" | tee {log}
          arcasHLA genotype {input.fq1} {input.fq2} -g A,B,C,DQA1,DQB1,DRB1 -t {threads} -v -o {params.path} 2>> {log}
        '''

## Makes a copy of the arcasHLA genotype json (output from rule arashla_genotype) from the arcashla directory into the merge directory.
## Currently does not have any further purpose in the pipeline
rule arcashla_relocate:
    input:
        json=rules.aracshla_genotype.output.json_genotype
    output:
        merge=paths.arcashla.merge
    benchmark:
        'benchmark/{sample}_arcashla_relocate.tab'
    log:
        'log/{sample}_arcashla_relocate.log'
    shell:
        '''
          echo "cp {input.json} {output.merge}" | tee {log}
          cp {input.json} {output.merge} 2>> {log}
        '''  
