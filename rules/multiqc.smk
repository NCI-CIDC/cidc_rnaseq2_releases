## Generates MultiQC report for STAR, Salmon, RSeQC, and Picard outputs
rule multiqc:
    input:
        star=rules.run_star.output,
      #  salmon=rules.salmon.output,
      #  bam_stat=rules.bam_stat.output,
      #  read_gc=rules.read_gc.output,
      #  tin_score=rules.tin_score.output,
      #  read_distribution=rules.read_distribution.output,
      #  gene_body_coverage=rules.gene_body_coverage.output,
      #  junction_saturation=rules.junction_saturation.output,
      #  collect_insert_size=rules.collect_insert_size.output
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
        dir=PREDIR+"/bam/*._STARgenome",
        dir2=PREDIR+"/bam/*._STARpass1",
        output=PREDIR+"/multiqc"
    shell:
        '''
          echo "multiqc {params.star} {params.salmon} {params.rseqc} --ignore {params.dir} {params.dir2} -o {params.output}" | tee {log}
          multiqc {params.star} {params.salmon} {params.rseqc} --ignore {params.dir} {params.dir2} -o {params.output} 2>> {log}
         
          ## Export rule env details
          conda env export --no-builds > info/multiqc.info 
        '''
