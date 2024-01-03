_microbiome_threads = 60

rule contamination_centrifuge_index:
    output:
        tar = paths.centrifuge.tar,
        db = paths.centrifuge.db
    log:
       to_log(paths.centrifuge.tar)
    message:
       "Building Centrifuge Index"
    benchmark:
       to_benchmark(paths.centrifuge.tar)
    threads: _microbiome_threads
    conda: "../envs/contamination.yaml"
    params:
       dest = Path(paths.centrifuge.tar).parent,
       URI = CFUG_REF
    shell:
       # '''curl -o {output.tar} https://genome-idx.s3.amazonaws.com/centrifuge/p_compressed%2Bh%2Bv.tar.gz;'''
        '''gsutil cp {params.URI} {params.dest} && '''
        '''tar -xvzf {output.tar} -C centrifuge'''

rule contamination_centrifuge:
     input:
        r1=expand(paths.rqual_filter.qfilter_fastq_paired, read=ENDS, paired=['P','U'])[0],
        r2=expand(paths.rqual_filter.qfilter_fastq_paired, read=ENDS, paired=['P','U'])[2],
        index=rules.contamination_centrifuge_index.output
     output:
        classification = paths.centrifuge.classification
     log:
        to_log(paths.centrifuge.classification)
     message:
        "Running Centrifuge on {wildcards.sample}"
     benchmark:
        to_benchmark(paths.centrifuge.classification)
     threads: _microbiome_threads
     conda: "../envs/contamination.yaml"
     params:
        tar_file=(rules.contamination_centrifuge_index.output[0]).replace('.tar.gz','')

	#if len(ENDS) == 2 else expa\
        #nd(paths.rqual_filter.qfilter_fastq_single, read=ENDS)[0]
     shell:
        """centrifuge -x {params.tar_file} -p {threads}  -q --host-taxids 9606 -1 {input.r1} -2 {input.r2} -S {output.classification}"""
