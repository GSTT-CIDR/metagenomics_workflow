rule delete_fastq:
    input:
        report = "/mnt/reports/{sample}/{sample}_{time}_hours_report.pdf",
        raw_fastq = "/mnt/results/{sample}/{time}_hours/files/{sample}_{time}_hours_merged.fastq"
    output:
        "/mnt/reports/{sample}/.{sample}_{time}_hours_report.pdf.tmp"
    run:
        shell("rm {input.raw_fastq}")
        shell("touch {output}")

rule compress_outputs:
    input:
        report = "/mnt/reports/{sample}/{sample}_{time}_hours_report.pdf",
        scrubbed = "/mnt/results/{sample}/{time}_hours/microbial/{sample}_{time}_hours_hg38_removed.fastq",
        centrifuge_raw = "/mnt/results/{sample}/{time}_hours/centrifuge/centrifuge_raw.tsv"
    output:
        compressed_scrubbed = "/mnt/results/{sample}/{time}_hours/microbial/{sample}_{time}_hours_hg38_removed.fastq.gz",
        compressed_centrifuge_raw = "/mnt/results/{sample}/{time}_hours/centrifuge/centrifuge_raw.tsv.gz"
    run:
        shell("pigz {input.scrubbed} ; pigz {input.centrifuge_raw}")
