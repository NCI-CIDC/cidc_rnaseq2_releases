## Set up directory structure based on dirs supplied in config
## Ignores non-zero exit status returned when any directories already exist
rule directory_setup:
    output:
        'progress/dirs.done'
    params:
        subdirs=SUBDIRS
    threads:1
    shell:
        '''
          mkdir {params.subdirs} -p 2> /dev/null
          touch {output}
        '''

## Download reference genome fa and supporting gtf annotation from ncbi ftp site
rule retrieve_reference_genome:
    input:
        rules.directory_setup.output
    output:
        fa=paths.genome.fa,
        gtf=paths.annot.gtf
    benchmark:
        'benchmark/retrieve_reference_genome.tab'
    log:
        'log/retrieve_reference_genome.log'
    params:
        fa_uri=GENOME_FA_URI,
        gtf_uri=GENOME_GTF_URI
    priority: 1000
    threads: 1
    shell:
        '''
          echo "downloading genome to map reads to GRCh38 or hg38..." | tee {log}
          gsutil cp {params.fa_uri} {output.fa}.gz && gunzip {output.fa}.gz
          
          echo "downloading supporting GTF annotations..." | tee -a {log}
          gsutil cp {params.gtf_uri} {output.gtf}.gz && gunzip {output.gtf}.gz
        '''

## Download built star_index files for the specified genome
## If using different genome, need to edit rule to call building index using STAR
rule build_star_index:
    input:
        rules.directory_setup.output,
        rules.retrieve_reference_genome.output.fa
    output:
        'progress/star_index_built.done'
    benchmark:
        'benchmark/build_star_index.tab'
    log:
        'log/build_star_index.log'
    conda:
        SOURCEDIR+"/../envs/star.yaml"
    params:
        star_uri=GENOME_STAR_URI
    priority: 1000
    threads: 1
    shell:
        '''
          echo "Downloading star_index files for mapping reads to GRCh38 or hg38..." | tee {log}
          gsutil -m cp -R {params.star_uri} genome
          touch {output} 
          
          ## export rule env details
          conda env export --no-builds > info/star.info
        '''

## Get genome chrom sizes for bw generation
rule genome_size:
    input:
        genome_fa=rules.retrieve_reference_genome.output.fa
    output:
        size=paths.genome.size,
        fai=paths.genome.fai
    benchmark:
        'benchmark/genome_size.tab'
    log:
        'log/genome_size.log'
    conda:
        SOURCEDIR+"/../envs/samtools.yaml"
    threads: 1
    shell:
        '''
          ## get genome chrom size
          echo "samtools faidx {input.genome_fa}" | tee {log}
          samtools faidx {input.genome_fa} 2>> {log}
          echo "cut -f1,2 {input.genome_fa}.fai > {output.size}" | tee -a {log}
          cut -f1,2 {input.genome_fa}.fai > {output.size} 2>> {log}
          
          ## export rule env details
          conda env export --no-builds > info/samtools.info
        '''

# Filter the hg38 genome index and convert from fai to bed format 
rule create_bed:
    input:
        paths.genome.fai
    output:
        paths.genome.bed
    benchmark:
        'benchmark/create_bed.tab'
    threads: 1
    shell:
        '''
          ## Remove the entries from chrM, chrUN, _random, chrEBV in the hg38 genome index and convert fai format to bed
          grep -v -E 'chrUn|_random|chrEBV|chrM' {input} | awk -F'\t' '{{ print $1,\"0\",$2 }}' > {output}
        '''  

## Retrieve hg38 blacklist from https://github.com/Boyle-Lab/Blacklist
rule retrieve_hg38_blacklist:
    output:
        paths.genome.blacklist
    benchmark:
        'benchmark/retrieve_hg38_blacklist.tab'
    params:
        blacklist_uri=GENOME_BLACKLIST_URI
    threads: 1
    shell:
        '''
          gsutil cp {params.blacklist_uri} {output}.gz
          gunzip {output}.gz
        '''

## Retrieve DHS regions list from dev GCP bucket. This might not be final location of the file.
## If file location changes, the shell directive needs to be updated.
rule retrieve_hg38_dhs:
    output:
        paths.genome.dhs
    benchmark:
        'benchmark/retrieve_hg38_dhs.tab'
    params:
        dhs_uri=GENOME_DHS_URI
    threads: 1
    shell:
        "gsutil cp {params.dhs_uri} {output}"

## Retrieve evolutionary bigwig file dev GCP bucket. This might not be final location of the file.
## If file location changes, the shell directive needs to be updated.
rule retrieve_conservation_bw:
    output:
        paths.annot.bw
    benchmark:
        to_benchmark(paths.annot.bw)
    params:
        dhs_uri=GENOME_CONSERVATION_URI
    threads: 1
    shell:
        "gsutil cp {params.dhs_uri} {output}"

## Retrieve hg38 RefSeq genes bed and hg38 housekeeping genes bed for the RSeQC module
rule retrieve_rseqc_beds:
    output:
       bed=paths.annot.refseq_bed,
       housekeeping_bed=paths.annot.housekeeping_bed
    benchmark:
       'benchmark/retrieve_rseqc_beds.tab'
    log:
       'log/retrieve_rseqc_beds.log'
    params:
       bed_uri=RSEQC_BED_URI,
       housekeeping_uri=RSEQC_HOUSEKEEPING_BED_URI
    shell:
       '''
       echo "gsutil cp {params.bed_uri} {output.bed}.gz && gunzip -k {output.bed}.gz" | tee {log}
       gsutil cp {params.bed_uri} {output.bed}.gz && gunzip -k {output.bed}.gz 2>> {log}

       echo "gsutil cp {params.housekeeping_uri} {output.housekeeping_bed}.gz && gunzip -k {output.housekeeping_bed}.gz" | tee -a {log}
       gsutil cp {params.housekeeping_uri} {output.housekeeping_bed}.gz && gunzip -k {output.housekeeping_bed}.gz 2>> {log}
       '''

## Retrieve Trinity Cancer Transcriptome Analysis Toolkit (CTAT) human hg38 library for STAR-Fusion 
rule retrieve_ctat_library:
    output:
        lib=directory(paths.annot.lib),
        tch='progress/ctat_lib_downloaded.done'
    benchmark:
        'benchmark/retrieve_ctat_library.tab'
    log:
        'log/retrieve_ctat_library.log'
    params:
        lib=FUSION_LIB_URI,
        predir=PREDIR
    shell:
        '''
          echo "gsutil -m cp -R {params.lib} annot && touch {output.tch}" | tee {log}
          gsutil -m cp -R {params.lib} annot && touch {output.tch} 2>> {log}
        '''

## Retrieve immune reference datasets
rule retrieve_immune_refs:
    output:
        response=paths.immune_response.msisensor2_models,
        repertoire=paths.immune_repertoire.bcrtcr
    benchmark:
        'benchmark/retrieve_immune_refs.tab'
    params:
        response_uri=IMMUNE_RESPONSE_URI,
        repertoire_uri=IMMUNE_REPERTOIRE_URI
    threads: 1
    shell:
        '''
        gsutil cp -r {params.response_uri} {output.response}
        gsutil cp {params.repertoire_uri} {output.repertoire}
        '''
