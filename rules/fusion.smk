## Identifies candidate fusion transcripts using STAR-Fusion
rule star_fusion:
    input:
        lib=rules.retrieve_ctat_library.output.lib,
        tch=rules.retrieve_ctat_library.output.tch,
        chim_junc=rules.run_star.output.chim_junc
    output:
        tsv=paths.fusion.tsv,
        tsv_abridged=paths.fusion.tsv_abridged,
        tsv_sample=paths.fusion.tsv_sample
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
          echo "STAR-Fusion --chimeric_junction {input.chim_junc} --genome_lib_dir {input.lib} --output_dir {params.predir}/fusion/{wildcards.sample} --CPU {threads} \
          && mv {params.predir}/fusion/{wildcards.sample}/star-fusion.fusion_predictions.tsv {output.tsv} \
          && mv {params.predir}/fusion/{wildcards.sample}/star-fusion.fusion_predictions.abridged.tsv {output.tsv_abridged} \
          awk 'BEGIN{{OFS="\t"}} NR==1{{$0=$0"\tsample"}} NR>1{{$0=$0"\t{wildcards.sample}"}} 1' {output.tsv_abridged} > {output.tsv_sample}" | tee {log}

          STAR-Fusion --chimeric_junction {input.chim_junc} --genome_lib_dir {input.lib} --output_dir {params.predir}/fusion/{wildcards.sample} --CPU {threads} \
          && mv {params.predir}/fusion/{wildcards.sample}/star-fusion.fusion_predictions.tsv {output.tsv} \
          && mv {params.predir}/fusion/{wildcards.sample}/star-fusion.fusion_predictions.abridged.tsv {output.tsv_abridged} \
          && awk 'BEGIN{{OFS="\t"}} NR==1{{$0=$0"\tsample"}} NR>1{{$0=$0"\t{wildcards.sample}"}} 1' {output.tsv_abridged} > {output.tsv_sample} 2>> {log} 

          ## Export rule env details
          conda env export --no-builds > info/star_fusion.info
        '''
