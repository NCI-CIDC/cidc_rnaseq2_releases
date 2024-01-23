## Identifies candidate fusion transcripts using STAR-Fusion
rule star_fusion:
    input:
        lib=rules.retrieve_ctat_library.output.lib,
        chim_junc=rules.run_star_output.chim_junc
    output:
        tsv=paths.fusion.tsv,
        tsv_abridged=paths.fusion.tsv_abridged
    benchmark:
        'benchmark/{sample}_star_fusion.tab'
    log:
        'log/{sample}_star_fusion.log'
    conda:
        SOURCEDIR+"/../envs/star_fusion.yaml"
    params:
        predir=PREDIR,
        
    threads: max(1,min(8,NCORES))
    shell:
        '''
          echo "STAR-Fusion --chimeric_junction {input.chim_junc} --genome_lib_dir {input.lib} --output_dir {params.predir}/fusion/{sample}_" | tee {log}
          STAR-Fusion --chimeric_junction {input.chim_junc} --genome_lib_dir {input.lib} --output_dir {params.predir}/fusion/{sample}_ --CPU {threads} 2>> {log}
        '''
