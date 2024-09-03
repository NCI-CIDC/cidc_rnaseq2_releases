## Generates MultiQC report for STAR, Salmon, RSeQC, and Picard outputs
rule multiqc:
    input:
        expand(rules.run_star.output, sample=SAMID),
        expand(rules.salmon.output, sample=SAMID),
        expand(rules.bam_stat.output, sample=SAMID),
        expand(rules.read_gc.output, sample=SAMID),
        expand(rules.tin_score.output, sample=SAMID),
        expand(rules.read_distribution.output, sample=SAMID),
        expand(rules.gene_body_coverage.output, sample=SAMID),
        expand(rules.junction_saturation.output, sample=SAMID),
        expand(rules.collect_insert_size.output, sample=SAMID)
    output:
        html=paths.multiqc.html
    benchmark:
        'benchmark/multiqc.tab'
    log:
        'log/multiqc.log'
    conda:
        SOURCEDIR+"/../envs/multiqc.yaml"
    params:
        star=PREDIR+"/bam",
        salmon=PREDIR+"/salmon",
        rseqc=PREDIR+"/rseqc",
        output=PREDIR+"/multiqc"
    shell:
        '''
          echo "multiqc {params.star} {params.salmon} {params.rseqc} -o {params.output}" | tee {log}
          multiqc {params.star} {params.salmon} {params.rseqc} -o {params.output} 2>> {log}
         
          ## Export rule env details
          conda env export --no-builds > info/multiqc.info 
        '''
