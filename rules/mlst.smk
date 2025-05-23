rule MLST:
    input:
        fastq = "/mnt/results/{sample}/{time}_hours/microbial/{sample}_{time}_hours_hg38_removed.fastq",
        read = "/mnt/results/{sample}/{time}_hours/centrifuge/read_assignments.tsv",
        targets = "/mnt/results/{sample}/{time}_hours/amr/centrifuge_top_hits.tsv"
    output:
        profiles = directory("/mnt/results/{sample}/{time}_hours/mlst/")
    shell:
        """
        mkdir {output.profiles}
        while read i; do
        echo "$i"
        out=$(echo "$i" | sed -e 's/ /_/g')
        path=$(python3 scripts/find_mlst_dir.py -t "$i" -d {config[mlst][directory]} -m {config[mlst][list]})
        if [ "$path" == None ]
        then
            echo "No MLST scheme" > {output.profiles}/${{out}}_mlst.tsv
        else
            python3 scripts/extract_reads.py -f {input.fastq} -c {input.read} -t "$i" | krocus $path - | tail -n 1 > {output.profiles}/${{out}}_mlst.tsv
        fi
        done < {input.targets}
        """
    



    