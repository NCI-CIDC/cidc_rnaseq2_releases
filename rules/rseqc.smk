## Summarizes mapping statistics of the aligned BAM
rule bam_stat:
    input:
        bam=rules.housekeeping_bam.output.bam if RSEQC=='housekeeping' else rules.downsample_bam.output.bam,
        bai=rules.housekeeping_bam.output.bai if RSEQC=='housekeeping' else rules.downsample_bam.output.bai
    output:
        txt=paths.rseqc.stat_txt
    benchmark:
        'benchmark/{sample}_bam_qc.tab'
    log:
        'log/{sample}_bam_qc.tab'
    conda:
        SOURCEDIR+"/../envs/rseqc.yaml"
    shell:
        '''
          echo "bam_stat.py -i {input.bam} > {output.txt}" | tee {log}
          bam_stat.py -i {input.bam} > {output.txt} 2>> {log}

          ## Export rule env details
          conda env export --no-builds > info/rseqc.info
        '''

## Calculates GC content distribution of reads with RSeQC 
rule read_gc:
    input:
        bam=rules.housekeeping_bam.output.bam if RSEQC=='housekeeping' else rules.downsample_bam.output.bam,
        bai=rules.housekeeping_bam.output.bai if RSEQC=='housekeeping' else rules.downsample_bam.output.bai
    output:
        r=paths.rseqc.gc_r,
        txt=paths.rseqc.gc_txt
    benchmark:
        'benchmark/{sample}_read_gc.tab'
    log:
        'log/{sample}_read_gc.log'
    conda:
        SOURCEDIR+"/../envs/rseqc.yaml"
    params:
        sample='{sample}'
    shell:
      '''
        echo "read_GC.py -i {input.bam} -o rseqc/read_gc/{params.sample}" | tee {log}
        read_GC.py -i {input.bam} -o rseqc/read_gc/{params.sample} 2>> {log}

        ## R script to get txt output info
        echo "out=as.vector(summary(gc));dta = data.frame('{params.sample}',out[1],out[2],out[3],out[4],out[5],out[6]);write.table(dta,file='{output.txt}',sep="\t",row.names=F,col.names=F,quote=F);" >> {output.r}
        sed -i "s/pdf/png/g" {output.r}
        sed -i 's/main=""/main="{params.sample}"/g' {output.r}
        Rscript --vanilla --quiet {output.r}
      '''

## Measures transcript integrity number (TIN) with RSeQC
rule tin_score:
    input:
        bam=rules.housekeeping_bam.output.bam if RSEQC=='housekeeping' else rules.downsample_bam.output.bam,
        bai=rules.housekeeping_bam.output.bai if RSEQC=='housekeeping' else rules.downsample_bam.output.bai,
        bed=rules.retrieve_rseqc_beds.output.housekeeping_bed if RSEQC=='housekeeping' else rules.retrieve_rseqc_beds.output.bed
    output:
        txt=paths.rseqc.ts_txt_hk if RSEQC=='housekeeping' else paths.rseqc.ts_txt,
        xls=paths.rseqc.ts_xls_hk if RSEQC=='housekeeping' else paths.rseqc.ts_xls
    benchmark:
        'benchmark/{sample}_tin_score.tab'
    log:
        'log/{sample}_tin_score.log'
    conda:
        SOURCEDIR+"/../envs/rseqc.yaml"
    params:
        predir=PREDIR,
        variable='housekeeping' if RSEQC=='housekeeping' else 'downsampling'
    shell:
        '''
          echo "tin.py --input={input.bam} --refgene={input.bed} --minCov=10 --sample-size=100" | tee {log}
          ## Default options used for minCov and sample-size
          tin.py --input={input.bam} --refgene={input.bed} --minCov=10 --sample-size=100 2>> {log}
          mv {params.predir}/{wildcards.sample}_{params.variable}.summary.txt {params.predir}/{wildcards.sample}_{params.variable}.tin.xls {params.predir}/rseqc/tin_score/
        ''' 

## Calculates how mapped reads were distributed over genome feature (CDS exon, 5’UTR exon, 3’ UTR exon, Intron, Intergenic regions) with RSeQC
rule read_distribution:
    input:
        bam=rules.housekeeping_bam.output.bam if RSEQC=='housekeeping' else rules.downsample_bam.output.bam,
        bai=rules.housekeeping_bam.output.bai if RSEQC=='housekeeping' else rules.downsample_bam.output.bai,
        bed=rules.retrieve_rseqc_beds.output.housekeeping_bed if RSEQC=='housekeeping' else rules.retrieve_rseqc_beds.output.bed
    output:
        txt=paths.rseqc.rd_txt
    benchmark:
        'benchmark/{sample}_read_distribution.tab'
    log:
        'log/{sample}_read_distribution.log'
    conda:
        SOURCEDIR+"/../envs/rseqc.yaml"
    shell:
        '''
          echo "read_distribution.py --input-file={input.bam} --refgene={input.bed} > {output.txt}" | tee {log}
          read_distribution.py --input-file={input.bam} --refgene={input.bed} > {output.txt} 2>> {log}
        '''

## Calculates the RNA-seq reads coverage over gene body with RSeQC. Requires input sorted bam and index
rule gene_body_coverage:
    input:
        bam=rules.housekeeping_bam.output.bam if RSEQC=='housekeeping' else rules.downsample_bam.output.bam,
        bai=rules.housekeeping_bam.output.bai if RSEQC=='housekeeping' else rules.downsample_bam.output.bai,
        bed=rules.retrieve_rseqc_beds.output.housekeeping_bed if RSEQC=='housekeeping' else rules.retrieve_rseqc_beds.output.bed
    output:
        txt=paths.rseqc.gbc_txt,
        r=paths.rseqc.gbc_r,
        png=paths.rseqc.gbc_png
    benchmark:
        'benchmark/{sample}_gene_body_coverage.tab'
    log:
        'log/{sample}_gene_body_coverage.log'
    conda:
        SOURCEDIR+"/../envs/rseqc.yaml"
    params:
        predir=PREDIR
    shell:
        '''
          echo "geneBody_coverage.py --input={input.bam} --refgene={input.bed} --format=png --out-prefix={params.predir}/rseqc/gene_body_coverage/{wildcards.sample}" | tee {log}
          geneBody_coverage.py --input={input.bam} --refgene={input.bed} --format=png --out-prefix={params.predir}/rseqc/gene_body_coverage/{wildcards.sample} 2>> {log}
        '''

## Checks junction saturation using RSeqQC 
rule junction_saturation:
    input:
        bam=rules.housekeeping_bam.output.bam if RSEQC=='housekeeping' else rules.downsample_bam.output.bam,
        bai=rules.housekeeping_bam.output.bai if RSEQC=='housekeeping' else rules.downsample_bam.output.bai,
        bed=rules.retrieve_rseqc_beds.output.housekeeping_bed if RSEQC=='housekeeping' else rules.retrieve_rseqc_beds.output.bed
    output:
        plot=paths.rseqc.js_plot,
        r=paths.rseqc.js_r
    benchmark:
        'benchmark/{sample}_junction_saturation.tab'
    log:
        'log/{sample}_junction_saturation.log'
    conda:
        SOURCEDIR+"/../envs/rseqc.yaml"
    params:
        predir=PREDIR
    shell:
        '''
          echo "junction_saturation.py --input-file={input.bam} --refgene={input.bed} --out-prefix={params.predir}/rseqc/junction_saturation/{wildcards.sample}" | tee {log}
          junction_saturation.py --input-file={input.bam} --refgene={input.bed} --out-prefix={params.predir}/rseqc/junction_saturation/{wildcards.sample} 2>> {log}
        '''

## Provides insert size distribution and read orientation of paired-end libraries using Picard
rule collect_insert_size:
    input:
        bam=rules.housekeeping_bam.output.bam if RSEQC=='housekeeping' else rules.downsample_bam.output.bam,
        bai=rules.housekeeping_bam.output.bai if RSEQC=='housekeeping' else rules.downsample_bam.output.bai,
        fa=paths.ref_files.fa
    output:
        txt=paths.rseqc.is_txt,
        pdf=paths.rseqc.is_pdf
    benchmark:
        'benchmark/{sample}_collect_insert_size.tab'
    log:
        'log/{sample}_collect_insert_size.log'
    conda:
        SOURCEDIR+"/../envs/picard.yaml"
    params:
        predir=PREDIR
    shell:
        '''
          echo "picard CollectInsertSizeMetrics I={input.bam} R={input.fa} M=0.5 O={output.txt} H={output.pdf}" | tee {log}
          picard CollectInsertSizeMetrics I={input.bam} R={input.fa} M=0.5 O={output.txt} H={output.pdf} 2>> {log}

          ## Export rule env details
          conda env export --no-builds > info/picard.info
        '''

## Moves the log.txt that is generated with rule gene_body_coverage (geneBody_coverage.py) from PREDIR to the rseqc/gene_body_coverage directory 
rule gene_body_coverage_log:
    input:
        expand(paths.rseqc.gbc_png, sample=SAMID)
    output:
        log=paths.rseqc.log
    benchmark:
        'benchmark/gene_body_coverage_log.tab'
    log:
        'log/gene_body_coverage_log.log'
    params:
        predir=PREDIR
    shell:
        '''
          echo "mv PREDIR/log.txt {output.log}" | tee {log}
          mv {params.predir}/log.txt {output.log} 2>> {log}
        '''
