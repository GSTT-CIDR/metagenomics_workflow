rule centrifuge_bacteria:
    input:
        micro = "/mnt/results/{sample}/{time}_hours/microbial/{sample}_{time}_hours_hg38_removed.fastq"
        # "/mnt/results/{sample}/{time}_hours/microbial/hg38_unmapped.fastq"
    output:
        raw = "/mnt/results/{sample}/{time}_hours/centrifuge/centrifuge_raw.tsv",
        report = temp("/mnt/results/{sample}/{time}_hours/centrifuge/centrifuge_report_raw.tsv")
    shell:
        """
        centrifuge -p 4 --mm -x {config[parameters][centrifuge][index][cmg]} -q {input.micro} -S {output.raw} \
        --report-file {output.report}
        """


rule parse_centrifuge:
    input:
        file = "/mnt/results/{sample}/{time}_hours/centrifuge/centrifuge_raw.tsv",
        fastq = "/mnt/results/{sample}/{time}_hours/microbial/{sample}_{time}_hours_hg38_removed.fastq" 
    output:
        report = "/mnt/results/{sample}/{time}_hours/centrifuge/centrifuge_report.tsv",
        read = "/mnt/results/{sample}/{time}_hours/centrifuge/read_assignments.tsv",
        failed = "/mnt/results/{sample}/{time}_hours/centrifuge/failed_reads.json",
        multi = "/mnt/results/{sample}/{time}_hours/centrifuge/multi_read.json"
    script:
        "../scripts/centrifuge_multi_match.py"

rule split_bacterial_viral:
    input:
        report = "/mnt/results/{sample}/{time}_hours/centrifuge/centrifuge_report.tsv",
    output:
        bacterial = "/mnt/results/{sample}/{time}_hours/centrifuge/bacterial_centrifuge_report.tsv",
        viral = "/mnt/results/{sample}/{time}_hours/centrifuge/viral_centrifuge_report.tsv"
    script:
        "../scripts/split_orders.py"
        





