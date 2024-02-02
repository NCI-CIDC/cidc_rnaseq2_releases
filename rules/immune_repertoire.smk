rule trust4_repertoire:
    input:
        bam=rules.run_star.output.bam,
        bai=rules.index_bam.output
    output:
        'trust4/{sample}_cdr3.out'
    benchmark:
        'benchmark/{sample}_trust4_immune_repertoire.tab'
    log:
        'log/{sample}.trust4_immune_repertoire.log'
    conda:
        SOURCEDIR+"/../envs/immune_repertoire.yaml"
    params:
        sample='{sample}',
        prefix = "trust4/{sample}",
        trust4_bcrtcr=paths.immune_repertoire.bcrtcr,
        genome=rules.retrieve_reference_genome.output.fa
    threads: max(1,min(8,NCORES))
    shell:
        '''
          run-trust4 -f {params.trust4_bcrtcr} --ref {params.genome} -b {input.bam} -t {threads} -o {params.prefix}
          rm -f {params.prefix}*fq
        '''

rule cdr3_process:
    input:
        bam=rules.trust4_repertoire.output,
    output:
        'trust4/{sample}_cdr3.out.processed.txt'
    benchmark:
        'benchmark/{sample}_cdr3_immune_repertoire.tab'
    log:
        'log/{sample}.cdr3_immune_repertoire.log'
    params:
        tmp='trust4/{sample}_cdr3.out.processed.tmp'
    shell:
        """
          perl trust-simplerep.pl {input.bam} > {params.tmp}
          sed -ig '1,1s/#count/count/g' {params.tmp}
          awk '{{print FILENAME}}' {params.tmp} | awk '{{print $3}}' FS='\t' | paste {params.tmp} - | awk -F '\t' 'NR==1{{$9="sample"}}1' OFS='\t'> {output}
          rm {params.tmp}
        """
