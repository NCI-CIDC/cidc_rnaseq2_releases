## Map reads to the reference genome using STAR and output sorted bam.
rule run_star:
    input:
        tch=rules.build_star_index.output,
        fa=rules.qualityfilter.output
    output:
        bam=paths.bam.bam,
        bam_transcriptome=paths.bam.bam_transcriptome,
        chim_junc=paths.bam.chimeric_junction
    benchmark:
        'benchmark/{sample}_run_star.tab'
    log:
        'log/{sample}_run_star.log'
    conda:
        SOURCEDIR+"/../envs/star.yaml"
    params:
        sample='{sample}',
        in_fa_str=expand(paths.rqual_filter.qfilter_fastq_paired, read=ENDS, paired=['P','U'])[0] + ' ' + expand(paths.rqual_filter.qfilter_fastq_paired, read=ENDS, paired=['P','U'])[2] if len(ENDS) == 2 else expand(paths.rqual_filter.qfilter_fastq_single, read=ENDS)[0]
    priority: 4
    threads: STAR
    shell:
        '''
          echo "STAR --runThreadN {threads} --genomeDir ref_files/star_index \
          --readFilesIn {params.in_fa_str} --readFilesCommand zcat \
          --outSAMstrandField intronMotif --outSAMunmapped Within --outSAMtype BAM SortedByCoordinate \
          --outReadsUnmapped None --chimSegmentMin 12 --chimJunctionOverhangMin 12 --chimOutJunctionFormat 1 \
          --alignSJDBoverhangMin 10 --alignMatesGapMax 100000 --alignIntronMax 100000 --alignSJstitchMismatchNmax 5 -1 5 5 \
          --chimMultimapScoreRange 3 --chimMultimapNmax 20 --chimNonchimScoreDropMin 10 \
          --peOverlapNbasesMin 12 --peOverlapMMp 0.1 \
          --twopassMode Basic --quantMode TranscriptomeSAM GeneCounts \
          --outSAMattrRGline ID:GRPundef SM:{wildcards.sample} \
          --outFileNamePrefix bam/{params.sample}. \
          --outStd BAM_SortedByCoordinate \
          --chimScoreJunctionNonGTAG -4 \
          --alignInsertionFlush Right \
          --alignSplicedMateMapLminOverLmate 0 \
          --alignSplicedMateMapLmin 30 > {output.bam}" | tee {log}

          STAR --runThreadN {threads} --genomeDir ref_files/star_index \
          --readFilesIn {params.in_fa_str} --readFilesCommand zcat \
          --outSAMstrandField intronMotif --outSAMunmapped Within --outSAMtype BAM SortedByCoordinate \
          --outReadsUnmapped None --chimSegmentMin 12 --chimJunctionOverhangMin 12 --chimOutJunctionFormat 1 \
          --alignSJDBoverhangMin 10 --alignMatesGapMax 100000 --alignIntronMax 100000 --alignSJstitchMismatchNmax 5 -1 5 5 \
          --chimMultimapScoreRange 3 --chimMultimapNmax 20 --chimNonchimScoreDropMin 10 \
          --peOverlapNbasesMin 12 --peOverlapMMp 0.1 \
          --twopassMode Basic --quantMode TranscriptomeSAM GeneCounts \
          --outSAMattrRGline ID:GRPundef SM:{wildcards.sample} \
          --outFileNamePrefix bam/{params.sample}. \
          --outStd BAM_SortedByCoordinate \
          --chimScoreJunctionNonGTAG -4 \
          --alignInsertionFlush Right \
          --alignSplicedMateMapLminOverLmate 0 \
          --alignSplicedMateMapLmin 30 > {output.bam} 2>> {log}

          ## Export rule env details
          conda env export --no-builds > info/star.info
        '''

## Index BAM
rule index_bam:
    input:
        bam=rules.run_star.output.bam
    output:
        paths.bam.index
    benchmark:
        'benchmark/{sample}_index_bam.tab'
    log:
        'log/{sample}_index_bam.log'
    conda:
        SOURCEDIR+"/../envs/samtools.yaml"
    params:
        sample='{sample}'
    priority: 5
    threads: max(1,min(8,NCORES))
    shell:
        '''
          echo "samtools index -@ {threads} {input.bam}" > {log}
          samtools index -@ {threads} {input.bam} 2>> {log}
        '''

## Generate stats for the raw aligned bam
rule align_bam_stats:
   input:
       bam=rules.run_star.output.bam,
       idx=rules.index_bam.output
   output:
       paths.bam.stats
   benchmark:
       'benchmark/{sample}_align_bam_stats.tab'
   log:
       'log/{sample}_align_bam_stats.log'
   conda:
       SOURCEDIR+"/../envs/samtools.yaml"
   threads: max(1,min(8,NCORES))
   shell:
       '''
         echo "samtools stats -@ {threads} {input.bam} | grep ^SN | cut -f 2- > {output}" | tee {log}
         samtools stats -@ {threads} {input.bam} | grep ^SN | cut -f 2- > {output} 2>> {log}
       '''

## Subsample the aligned bam for use in the RSeQC module; the output bam is sorted
rule downsample_bam:
    input:
        stats=rules.align_bam_stats.output,
        bam=rules.run_star.output.bam,
        idx=rules.index_bam.output
    output:
        seq=paths.bam.size,
        bam=paths.bam.downsampled_bam,
        bai=paths.bam.downsampled_bai
    benchmark:
        'benchmark/{sample}_downsample_bam.tab'
    log:
        'log/{sample}_downsample_bam.log'
    conda:
        SOURCEDIR+"/../envs/samtools.yaml"
    params:
        srcdir=SOURCEDIR
    threads: max(1,min(8,NCORES))
    shell:
        '''
          echo "bash {params.srcdir}/shell/downsample_bam.sh {input.stats} {output.seq} {input.bam} {output.bam} {output.bai} {threads}" | tee {log}
          bash {params.srcdir}/shell/downsample_bam.sh {input.stats} {output.seq} {input.bam} {output.bam} {output.bai} {threads} 2>> {log}
        '''

## Generate bam with only the housekeeping genes from the subsampled bam for use in the RSeQC module
rule housekeeping_bam:
    input:
        bam=rules.downsample_bam.output.bam,
        bai=rules.downsample_bam.output.bai,
        stats=rules.downsample_bam.output.seq,
        bed=rules.retrieve_rseqc_beds.output.housekeeping_bed
    output:
        bam=paths.bam.housekeeping_bam,
        bai=paths.bam.housekeeping_bai
    benchmark:
        'benchmark/{sample}_housekeeping_bam.tab'
    log:
        'log/{sample}_housekeeping_bam.log'
    conda:
        SOURCEDIR+"/../envs/bedtools.yaml"
    threads: max(1,min(8,NCORES))
    shell:
        '''
          echo "bedtools intersect -a {input.bam} -b {input.bed} > {output.bam} && samtools index -@ {threads} {output.bam} > {output.bai}" | tee {log}
          bedtools intersect -a {input.bam} -b {input.bed} > {output.bam} && samtools index -@ {threads} {output.bam} > {output.bai} 2>> {log}

          ## Export rule env details
          conda env export --no-builds > info/bedtools.info
        '''
