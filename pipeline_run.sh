#!/bin/bash
eval "$(conda shell.bash hook)"
conda activate cmg

read -p "Please enter sample table name, leave blank for default: [sample_table.tsv]: " sample
sample=${sample:-sample_table.tsv}
echo "Using sample table: $sample"
read -p "How many samples to run: " num
echo "Running with $num cores"
echo "Do you want to run pipeline or perform dry run (to check parameters)?"
select yn in "Run" "Check"; do
	case $yn in
		Check ) for t in {2}; do snakemake --cores $num -k --config time=$t samples=$sample -n; done; break;;
		Run ) for t in {0.5,1,2,16,24}; do snakemake --cores $num -k --config time=$t samples=$sample --latency-wait 15; done; break;;
	esac
done

#snakemake --cores 9 -k --config time=0.5 -n
#snakemake --cores 9 -k --config time=2 -n
#snakemake --cores 9 -k --config time=16 -n


conda deactivate

