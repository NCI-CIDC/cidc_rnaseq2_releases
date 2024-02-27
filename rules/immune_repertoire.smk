## Run TRUST4 to provide an overview of the immune repertoire of tumors including CDR3 sequence length and the frequency of various V genes, J genes, and VJ pairs for different chains in the TCR and BCR
rule trust4:
    input:
        bam=rules.run_star.output.bam,
        bai=rules.index_bam.output,
        bcrtcr=rules.retrieve_immune_refs.output.bcrtcr,
        imgt=rules.retrieve_immune_refs.output.imgt
    output:
        cdr3=paths.trust4.cdr3
    benchmark:
        'benchmark/{sample}_trust4.tab'
    log:
        'log/{sample}_trust4.log'
    conda:
        SOURCEDIR+"/../envs/trust4.yaml"
    params:
        prefix = "trust4/{sample}/{sample}"
    threads: max(1,min(8,NCORES))
    shell:
        '''
          echo "run-trust4 -f {input.bcrtcr} --ref {input.imgt} -b {input.bam} -t {threads} -o {params.prefix} --abnormalUnmapFlag \
          && rm -f {params.prefix}*fq" | tee {log}

          run-trust4 -f {input.bcrtcr} --ref {input.imgt} -b {input.bam} -t {threads} -o {params.prefix} --abnormalUnmapFlag \
          && rm -f {params.prefix}*fq 2>> {log}

          ## Export rule env details
          conda env export --no-builds > info/trust4.info
        '''

## Run TRUST4 simplerep.pl script. This rule was from the old pipeline; however, it might not have any functional use currently.
rule cdr3_preprocess:
    input:
        cdr3=rules.trust4.output.cdr3
    output:
        txt=paths.trust4.txt,
        tmp=temp('trust4/{sample}/{sample}_cdr3.out.processed.tmp'),
        tmpg=temp('trust4/{sample}/{sample}_cdr3.out.processed.tmpg')
    benchmark:
        'benchmark/{sample}_cdr3_preprocess.tab'
    log:
        'log/{sample}_cdr3_preprocess.log'
    conda:
        SOURCEDIR+"/../envs/trust4.yaml"
    shell:
        '''
          ## Identifies path to the simplerep.pl script from TRUST4 in the current Conda environment
          script_path=$(conda info --envs | grep -E '\*' | awk '{{print $NF}}')

          echo "perl ${{script_path}}/bin/trust-simplerep.pl {input.cdr3} > {output.tmp} \
          && sed -ig '1,1s/#count/count/g' {output.tmp} 2>> {log} \
          && awk '{{print FILENAME}}' {output.tmp} | awk '{{print \$3}}' FS='\t' | paste {output.tmp} - | awk -F '\t' 'NR==1{{\$9="sample"}}1' OFS='\t'> {output.txt}" | tee {log}

          perl ${{script_path}}/bin/trust-simplerep.pl {input.cdr3} > {output.tmp} \
          && sed -ig '1,1s/#count/count/g' {output.tmp} \
          && awk '{{print FILENAME}}' {output.tmp} | awk '{{print $3}}' FS='\t' | paste {output.tmp} - | awk -F '\t' 'NR==1{{$9="sample"}}1' OFS='\t'> {output.txt} 2>> {log}
        '''
