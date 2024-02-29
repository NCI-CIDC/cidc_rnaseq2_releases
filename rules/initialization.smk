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
        fa=paths.ref_files.fa,
        gtf=paths.ref_files.gtf
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
          echo "gsutil cp {params.fa_uri} {output.fa}.gz && gunzip {output.fa}.gz" | tee {log}
          gsutil cp {params.fa_uri} {output.fa}.gz && gunzip {output.fa}.gz
          
          echo "gsutil cp {params.gtf_uri} {output.gtf}.gz && gunzip {output.gtf}.gz" | tee -a {log}
          gsutil cp {params.gtf_uri} {output.gtf}.gz && gunzip {output.gtf}.gz
        '''

## Download built star_index files for the specified genome
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
    params:
        star_uri=GENOME_STAR_URI
    priority: 1000
    threads: 1
    shell:
        '''
          echo "gsutil -m cp -R {params.star_uri} ref_files" | tee {log}
          gsutil -m cp -R {params.star_uri} ref_files
          touch {output} 
        '''

## Retrieve Gencode v37 comprehesive genes bed and Gencode v37 housekeeping genes bed for the RSeQC module
rule retrieve_rseqc_beds:
    output:
       bed=paths.ref_files.comp_bed,
       housekeeping_bed=paths.ref_files.housekeeping_bed
    benchmark:
       'benchmark/retrieve_rseqc_beds.tab'
    log:
       'log/retrieve_rseqc_beds.log'
    params:
       bed_uri=RSEQC_BED_URI,
       housekeeping_uri=RSEQC_HOUSEKEEPING_BED_URI
    shell:
       '''
       echo "gsutil cp {params.bed_uri} {output.bed}.gz && gunzip {output.bed}.gz" | tee {log}
       gsutil cp {params.bed_uri} {output.bed}.gz && gunzip {output.bed}.gz 2>> {log}

       echo "gsutil cp {params.housekeeping_uri} {output.housekeeping_bed}.gz && gunzip {output.housekeeping_bed}.gz" | tee -a {log}
       gsutil cp {params.housekeeping_uri} {output.housekeeping_bed}.gz && gunzip {output.housekeeping_bed}.gz 2>> {log}
       '''

## Retrieve Trinity Cancer Transcriptome Analysis Toolkit (CTAT) human hg38 library for STAR-Fusion 
rule retrieve_ctat_library:
    output:
        lib=directory(paths.ref_files.lib),
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
          echo "gsutil -m cp -R {params.lib} ref_files && touch {output.tch}" | tee {log}
          gsutil -m cp -R {params.lib} ref_files && touch {output.tch} 2>> {log}
        '''

## Retrieve reference datasets for immune repsonse (MSIsensor2) and immume repertoire (TRUST4) modules
rule retrieve_immune_refs:
    output:
        models=directory(paths.ref_files.models),
        tch='progress/msisensor2_models_downloaded.done',
        bcrtcr=paths.ref_files.bcrtcr,
        imgt=paths.ref_files.imgt
    benchmark:
        'benchmark/retrieve_immune_refs.tab'
    log:
        'log/retrieve_immune_refs.log'
    params:
        response_uri=IMMUNE_RESPONSE_URI,
        repertoire_uri=IMMUNE_REPERTOIRE_URI,
        repertoire_imgt_uri=IMMUNE_REPERTOIRE_IMGT_URI
    shell:
        '''
          echo "gsutil -m cp -R {params.response_uri} ref_files && touch {output.tch}" | tee {log}
          gsutil -m cp -R {params.response_uri} ref_files && touch {output.tch} 2>> {log}
           
          echo "gsutil cp {params.repertoire_uri} {output.bcrtcr} && gsutil cp {params.repertoire_imgt_uri} {output.imgt}" | tee -a {log}
          gsutil cp {params.repertoire_uri} {output.bcrtcr} && gsutil cp {params.repertoire_imgt_uri} {output.imgt} 2>> {log}
        '''

## Retrieve Centrifuge index (bacteria, archaea, viruses, and human) for the microbiome module
rule retrieve_centrifuge_idx:
    output:
        tar=temp(paths.ref_files.tar),
        idx=paths.ref_files.idx,
        tch='progress/centrifuge_idx_downloaded.done'
    benchmark:
        'benchmark/retrieve_centrifuge_idx.tab'
    log:
        'log/retrieve_centrifuge_idx.log'
    params:
        cfug_uri=CFUG_REF
    shell:
        '''
          echo "gsutil cp {params.cfug_uri} {output.tar} && tar -xvzf {output.tar} -C ref_files \
          && touch {output.tch}" | tee {log}

          gsutil cp {params.cfug_uri} {output.tar} && tar -xvzf {output.tar} -C ref_files \
          && touch {output.tch} 2>> {log}
        '''

## Retrieve Gencode transcripts fasta for Salmon in the quantification module
rule retrieve_transcripts_fa:
    output:
        fa=paths.ref_files.transcripts
    benchmark:
        'benchmark/retrieve_transcripts_fa.tab'
    log:
        'log/retrieve_transcripts_fa.log'
    params:
        fa_uri=TRANSCRIPTS_FA_URI
    shell:
        '''
          echo "gsutil cp {params.fa_uri} {output.fa}" | tee {log}
          gsutil cp {params.fa_uri} {output.fa} 2>> {log}
        '''

## Clone the IPD-IMGT/HLA repository and set up the database for use in the neoantigen module (arcasHLA).
## Due to the increasing size of the hla.dat file, the repository now requires the use of the Git LFS 
## tools (https://git-lfs.github.com) to handle files over 100MB in size.
rule retrieve_imgthla_db:
    output:
        tch='progress/imgthla_db_downloaded.done'
    benchmark:
        'benchmark/retrieve_imgthla_db.tab'
    log:
        'log/retrieve_imgthla_db.log'
    conda:
        SOURCEDIR+"/../envs/arcashla.yaml"
    shell:
        '''
          echo "git lfs install && arcasHLA reference --commit df6ba6f --verbose && touch {output.tch}" | tee {log}
          git lfs install && arcasHLA reference --commit df6ba6f --verbose && touch {output.tch} 2>> {log}
        '''
