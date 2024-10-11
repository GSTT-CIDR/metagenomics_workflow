rule nanostat:
    input:
        "/mnt/results/{sample}/{time}_hours/microbial/{sample}_{time}_hours_hg38_removed.fastq"
    output:
        "/mnt/results/{sample}/{time}_hours/qc/nanostat_summary.txt",
 
    shell:
        "NanoStat --fastq {input} -n {output}"

rule sample_table:
    output:
        "/mnt/results/{sample}/{time}_hours/qc/sample_table.tsv"
    shell:
        "cp {config[samples]} > {output}"

#rule get_versions:
#    output:
#        "/mnt/results/{sample}/{time}_hours/qc/versions.txt"
#    shell:
#    """
    
    


