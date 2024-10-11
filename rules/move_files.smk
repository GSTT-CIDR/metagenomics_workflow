rule move:
    input:
        seq_dir = lambda wildcards: sample_table.path[wildcards.sample]
    output:
        analysis = "/mnt/results/{sample}/{time}_hours/files/{sample}_{time}_hours_merged.fastq"
    log:
        "/mnt/results/{sample}/{time}_hours/log/file_move.log"
    script:
        "../scripts/transfer_reads.py"

rule copy_sample_sheet:
    input:
        sample_sheet = sample_sheet_path
    output:
        sample_sheet = "/mnt/results/{sample}/{time}_hours/files/sample_sheet.tsv"
    shell:
        """
        # Copy sample_sheet to the new directory
        cp {input.sample_sheet} {output.sample_sheet}
        """
