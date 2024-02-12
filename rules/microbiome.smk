## Run Centrifuge to classify sequences by assigning them to taxonomic categories 
rule centrifuge:
    input:
        r1=expand(paths.rqual_filter.qfilter_fastq_paired, read=ENDS, paired=['P','U'])[0],
        r2=expand(paths.rqual_filter.qfilter_fastq_paired, read=ENDS, paired=['P','U'])[2],
        idx=rules.retrieve_centrifuge_idx.output.idx,
        tch=rules.retrieve_centrifuge_idx.output.tch
    output:
        gz=paths.centrifuge.gz,
        tsv=paths.centrifuge.tsv,
        tsv_sample=paths.centrifuge.tsv_sample 
    benchmark:
        'benchmark/{sample}_centrifuge.tab'
    log:
        'log/{sample}_centrifuge.log'
    conda: 
        SOURCEDIR+"/../envs/centrifuge.yaml"
    params:
        idx=PREDIR+"/genome/p_compressed+h+v",
        txt=paths.centrifuge.txt
    threads: max(1,min(8,NCORES))
    shell:
        '''
          echo "centrifuge -x {params.idx} -p {threads} --host-taxids 9606 -1 {input.r1} -2 {input.r2} -S {params.txt} --report-file {output.tsv} \
          && gzip {params.txt} \
          && awk -v sampleID="{wildcards.sample}" 'BEGIN {{OFS="\t"}} {{if (NR == 1) print "sample", $0; else print sampleID, $0}}' {output.tsv} > {output.tsv_sample}" | tee {log}

          centrifuge -x {params.idx} -p {threads} --host-taxids 9606 -1 {input.r1} -2 {input.r2} -S {params.txt} --report-file {output.tsv} \
          && gzip {params.txt} \
          && awk -v sampleID="{wildcards.sample}" 'BEGIN {{OFS="\t"}} {{if (NR == 1) print "sample", $0; else print sampleID, $0}}' {output.tsv} > {output.tsv_sample} 2>> {log}
        '''
