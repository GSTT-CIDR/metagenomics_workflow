rule human_reads:
    input:
        "/mnt/results/{sample}/{time}_hours/files/{sample}_{time}_hours_merged.fastq"
    output:
        temp("/mnt/results/{sample}/{time}_hours/host/{sample}_{time}_hours_hg38_mapped.sam")
    shell:
        "minimap2 -x map-ont --secondary no -a {config[parameters][hg38][index]} {input} > {output}"

rule micro_fastq:
    input:
        "/mnt/results/{sample}/{time}_hours/host/{sample}_{time}_hours_hg38_mapped.sam"
    output:
        "/mnt/results/{sample}/{time}_hours/microbial/{sample}_{time}_hours_hg38_removed.fastq"
        # "/mnt/results/{sample}/{time}_hours/microbial/hg38_unmapped.fastq"
    shell:
        "samtools fastq -f 4 {input} > {output}"

rule map_stats:
    input:
        "/mnt/results/{sample}/{time}_hours/host/{sample}_{time}_hours_hg38_mapped.sam"
    output:
        "/mnt/results/{sample}/{time}_hours/host/{sample}_{time}_hours_map_stats.txt"
    shell:
        "samtools flagstat {input} > {output}"

