rule abricate_amr:
    input:
        fq = "/mnt/results/{sample}/{time}_hours/microbial/{sample}_{time}_hours_hg38_removed.fastq"
    output:
        amr = "/mnt/results/{sample}/{time}_hours/amr/amr_results.tsv"
    shell:
        "abricate --quiet --mincov 90 --db card {input.fq} > {output.amr}"


rule top_centrifuge:
    input:
        "/mnt/results/{sample}/{time}_hours/centrifuge/centrifuge_report.tsv"
    output:
        "/mnt/results/{sample}/{time}_hours/amr/centrifuge_top_hits.tsv"
    script:
        "../scripts/scagaire_targets.py"


rule scagaire:
    input:
        amr_res = "/mnt/results/{sample}/{time}_hours/amr/amr_results.tsv",
        top_hit = "/mnt/results/{sample}/{time}_hours/amr/centrifuge_top_hits.tsv"
    output:
        summary = "/mnt/results/{sample}/{time}_hours/amr/scagaire_gene_summary.tsv",
        report = "/mnt/results/{sample}/{time}_hours/amr/scagaire_report.tsv"
    shell:
        """
        t=$(cat {input.top_hit} | paste -sd "," -)
        scagaire "$t" {input.amr_res} -n card -s {output.summary} -o {output.report}
        if [ ! -f "{output.summary}" ]; then
            echo "No species in database" > {output.summary}
            echo "No Report" > {output.report}
        fi
        """


rule abricate_virulence:
    input:
        fq = "/mnt/results/{sample}/{time}_hours/microbial/{sample}_{time}_hours_hg38_removed.fastq"
    output:
        vdf = "/mnt/results/{sample}/{time}_hours/amr/virulence_factor_raw.tsv"
    shell:
        "abricate --quiet --mincov 90 --db vfdb {input.fq} > {output.vdf}"


rule parse_virulence_results:
    input:
        vdf = "/mnt/results/{sample}/{time}_hours/amr/virulence_factor_raw.tsv"
    output:
        res = "/mnt/results/{sample}/{time}_hours/amr/virulence_factor_summary.tsv"
    script:
        "../scripts/virulence_summary.py"



