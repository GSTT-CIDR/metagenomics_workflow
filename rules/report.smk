rule compile_report:
    input:
        centrifuge = "/mnt/results/{sample}/{time}_hours/centrifuge/bacterial_centrifuge_report.tsv",
        centrifuge_raw = "/mnt/results/{sample}/{time}_hours/centrifuge/centrifuge_raw.tsv",
        viral = "/mnt/results/{sample}/{time}_hours/centrifuge/viral_centrifuge_report.tsv",
        amr_summary = "/mnt/results/{sample}/{time}_hours/amr/scagaire_gene_summary.tsv",
        amr_report ="/mnt/results/{sample}/{time}_hours/amr/scagaire_report.tsv",
        vf = "/mnt/results/{sample}/{time}_hours/amr/virulence_factor_summary.tsv",
        qc="/mnt/results/{sample}/{time}_hours/qc/nanostat_summary.txt",
        stats = "/mnt/results/{sample}/{time}_hours/host/{sample}_{time}_hours_map_stats.txt"
    output:
        "/mnt/reports/{sample}/{sample}_{time}_hours_report.pdf"

    script:
        "../scripts/generate_report.py"

