## Get input
rule getfile:
    input:
        rules.directory_setup.output
    output:
        expand(paths.input.input_fastq, read=ENDS) if not BAM_INPUT else paths.input.input_bam
    benchmark:
        'benchmark/{sample}_getfile.tab'
    log:
        'log/{sample}_getfile.log'
    conda:
        SOURCEDIR+"/../envs/getfile.yaml"
    params:
        sample='{sample}',
        srcdir=SOURCEDIR,
        ends=','.join(ENDS),
        samid=','.join(SAMID),
        fastq1=','.join(FASTQ_1),
        fastq2=','.join(FASTQ_2),
        bam=','.join(BAM),
        cloud=CLOUD,
        output_joined=','.join(expand(paths.input.input_fastq, read=ENDS)) if not BAM_INPUT else paths.input.input_bam
    priority: 2
    threads: max(1,min(8,NCORES))
    shell:
        '''
          echo python3 {params.srcdir}/python/init-getfile.py --src {params.srcdir} --ends {params.ends} \
          --samid {params.samid} --fastq1 {params.fastq1} --fastq2 {params.fastq2} --bam {params.bam} \
          --cloud {params.cloud} --fastqdump parallel-fastq-dump \
          --sample {params.sample} --output {params.output_joined} \
          > {log}

          python3 {params.srcdir}/python/init-getfile.py --src {params.srcdir} --ends {params.ends} \
          --samid {params.samid} --fastq1 {params.fastq1} --fastq2 {params.fastq2} --bam {params.bam} \
          --cloud {params.cloud} --fastqdump parallel-fastq-dump \
          --sample {params.sample} --output {params.output_joined} \
          2>> {log}

          ## export rule env details
          conda env export --no-builds > info/getfile.info
        '''

## TODO rule for converting bam to PE fastq
rule bam2fastq:
    input:
        bam=paths.input.input_bam
    output:
        fq1=expand(paths.input.input_fastq, read=['1']),
        fq2=expand(paths.input.input_fastq, read=['2']),
    benchmark:
        'benchmark/{sample}_bam2fastq.tab'
    log:
        'log/{sample}_bam2fastq.log'
    conda:
        SOURCEDIR+"/../envs/samtools.yaml"
    priority: 3
    threads: max(1,min(8,NCORES))
    shell:
        '''
          samtools collate -@ {threads} -u -O {input.bam} | \
          samtools fastq -@ {threads} -1 {output.fq1} -2 {output.fq2} -0 /dev/null -s /dev/null -n
        '''

### Optionally trim adapters with cutadapt
rule trimadapters:
    input:
        fa=expand(paths.input.input_fastq, read=ENDS)
    output:
        [x + TRIM_ADAPTERS_OUTPUT for x in expand(paths.cutadapt.cutadapt_fastq, read=ENDS)]
    benchmark:
        'benchmark/{sample}_trimadapters.tab'
    log:
        'log/{sample}_trimadapters.log'
    conda:
         SOURCEDIR+"/../envs/trimadpaters.yaml"
    params:
        sample='{sample}',
        srcdir=SOURCEDIR,
        samid=','.join(SAMID),
        trim_fp=TRIM_FP,
        trim_tp=TRIM_TP,
        fp_adapters=','.join(FP_ADAPTERS),
        tp_adapters=','.join(TP_ADAPTERS),
        input_joined=','.join(expand(paths.input.input_fastq, read=ENDS)),
        output_joined=','.join([x + TRIM_ADAPTERS_OUTPUT for x in expand(paths.cutadapt.cutadapt_fastq, read=ENDS)])
    priority: 3
    threads: max(1,min(8,NCORES))
    shell:
      '''
        echo python3 {params.srcdir}/python/init-trimadapters.py --src {params.srcdir} --trim-fp {params.trim_fp} --trim-tp {params.trim_tp}\
        --samid {params.samid} --fp-adapters {params.fp_adapters} --tp-adapters {params.tp_adapters} \
        --cutadapt cutadapt --sample {params.sample} --input {params.input_joined} --output {params.output_joined} --threads {threads}\
        > {log}

        python3 {params.srcdir}/python/init-trimadapters.py --src {params.srcdir} --trim-fp {params.trim_fp} --trim-tp {params.trim_tp}\
        --samid {params.samid} --fp-adapters {params.fp_adapters} --tp-adapters {params.tp_adapters} \
        --cutadapt cutadapt --sample {params.sample} --input {params.input_joined} --output {params.output_joined} --threads {threads} --log {log}\
        2>> {log}

        ## export rule env details
        conda env export --no-builds > info/trimadapters.info
      '''

## Optionally quality-trim reads with Trimmomatic
rule qualityfilter:
    input:
        rules.trimadapters.output if TRIM_FP or TRIM_TP else expand(paths.input.input_fastq, read=ENDS)
    output:
        expand(paths.rqual_filter.qfilter_fastq_paired, read=ENDS, paired=['P','U']) if len(ENDS)==2 else expand(paths.rqual_filter.qfilter_fastq_single, read=ENDS)
    benchmark:
        'benchmark/{sample}_qualityfilter.tab'
    log:
        'log/{sample}_qualityfilter.log'
    conda:
        SOURCEDIR+"/../envs/qualityfilter.yaml"
    params:
        sample='{sample}',
        srcdir=SOURCEDIR,
        ends=','.join(ENDS),
        qual_cutoff=QUAL_CUTOFF,
        input_joined=','.join([x + TRIM_ADAPTERS_OUTPUT for x in expand(paths.cutadapt.cutadapt_fastq, read=ENDS)] if TRIM_FP or TRIM_TP else expand(paths.input.input_fastq, read=ENDS)),
        output_joined=','.join(expand(paths.rqual_filter.qfilter_fastq_paired, read=ENDS, paired=['P','U']) if len(ENDS)==2 else expand(paths.rqual_filter.qfilter_fastq_single, read=ENDS))
    priority: 3
    threads: max(1,min(8,NCORES))
    shell:
      '''
        echo python3 {params.srcdir}/python/init-qualityfilter.py --src {params.srcdir} --ends {params.ends} --threads {threads} \
        --trimmomatic trimmomatic --qual-cutoff {params.qual_cutoff} --input {params.input_joined} --output {params.output_joined} --log {log}\
        > {log}

        python3 {params.srcdir}/python/init-qualityfilter.py --src {params.srcdir} --ends {params.ends} --threads {threads} \
        --trimmomatic trimmomatic --qual-cutoff {params.qual_cutoff} --input {params.input_joined} --output {params.output_joined} --log {log} \
        2>> {log}

        ## export rule env details
        conda env export --no-builds > info/qualityfilter.info
      '''
