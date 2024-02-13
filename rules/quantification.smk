_salmon_threads = 16

rule quantification_gffread:
    input:
        gtf = paths.annot.gtf,
        fna = paths.genome.fa
    output:
        target_fa = paths.gffread.target_fa
    log:
       to_log(paths.gffread.target_fa)
    message:
       "Building target for salmon"
    benchmark:
       to_benchmark(paths.gffread.target_fa)
    threads: _salmon_threads
    conda: "../envs/quantification.yml"
    shell:
        '''gffread -w {output.gffread} -g {input.fna} {input.gtf}'''

rule quantification_salmon:
     input:
         bam=paths.bam.bam,
         target_fa=rules.quantification_gffread.output
     output:
        quants = paths.salmon.quants
     log:
        to_log(paths.salmon.quants)
     message:
        "Running Salmon on {sample}"
     benchmark:
        to_benchmark(paths.salmon.quants)
     threads: _salmon_threads
     conda: "../envs/quantification.yml"
     params:
        folder = Path(paths.salmon.quants).parent,
        sample = "{sample}"
     shell:
        "salmon quant -a {input.bam} -p {threads} -o {params.folder}/{params.sample} -t {input.target_fa} -l A && "
        "mv {params.folder}/{params.sample}/quants.sf {output}"
