rule human_reads:
    input:
        raw_fastq = "/mnt/results/{sample}/{time}_hours/files/{sample}_{time}_hours_merged.fastq",
        centrifuge_raw = "/mnt/results/{sample}/{time}_hours/centrifuge/centrifuge_raw.tsv"
    output:
        all_read_ids = "/mnt/results/{sample}/{time}_hours/host/{sample}_{time}_hours_all_read_ids.txt",
        human_read_ids = "/mnt/results/{sample}/{time}_hours/host/{sample}_{time}_hours_human_read_ids.txt"

    shell:
        # New function replacing the hg38 mapping
        """seqkit seq -ni {input.raw_fastq} > {output.all_read_ids} ; tail -n +2 {input.centrifuge_raw} | sort -k4,4rn -k1,1 | awk '!seen[$1]++' | awk '$3 == 9606  {{print $1}}'  > {output.human_read_ids}"""


rule micro_fastq:
    input:
        read_ids = "/mnt/results/{sample}/{time}_hours/host/{sample}_{time}_hours_human_read_ids.txt",
        raw_fastq = "/mnt/results/{sample}/{time}_hours/files/{sample}_{time}_hours_merged.fastq"
    output:
        hg38_removed = "/mnt/results/{sample}/{time}_hours/microbial/{sample}_{time}_hours_hg38_removed.fastq"
    shell:
        """seqkit grep -vf {input.read_ids} {input.raw_fastq} -o {output.hg38_removed}"""


rule map_stats:
    input:
        all_read_ids = "/mnt/results/{sample}/{time}_hours/host/{sample}_{time}_hours_all_read_ids.txt",
        human_read_ids = "/mnt/results/{sample}/{time}_hours/host/{sample}_{time}_hours_human_read_ids.txt"
    output:
        combined_stats = "/mnt/results/{sample}/{time}_hours/host/{sample}_{time}_hours_map_stats.txt"

    shell:
        """
        total_reads=$(cat {input.all_read_ids} | wc -l | xargs)
        human_reads=$(cat {input.human_read_ids} | wc -l | xargs)
        human_percentage=$(echo "${{human_reads}} * 100 / ${{total_reads}}" | bc)
        # Print the total reads
        echo $total_reads > {output.combined_stats}
        # Print the human reads count and percentage
        echo "${{human_reads}}/${{human_percentage}}%" >> {output.combined_stats}

        """

