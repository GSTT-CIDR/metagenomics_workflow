import pandas as pd
import glob
import time
import os
import sys
configfile: "/mnt/configs/metagenomics_config.yaml"


def find_path(exp, sample, barcode, data_dir):
    barcode_dict = {"1": "01", "2":"02", "3":"03", "4":"04", "5":"05", "6":"06", "7":"07", "8":"08", "9":"09", "10":"10", "11":"11", "12a":"12a", "12":"12a"}
    path_str = f"{str(data_dir)}/{exp}/{sample}/**/fastq_pass/barcode{barcode_dict[str(barcode)]}/"
    path = glob.glob(path_str, recursive=True)
    path = sorted(path, reverse=True)
    if len(path) == 0:
        sys.exit(f"Error: path error with sample {exp} barcode {barcode}. Check variables for experiment, sample or barcode")
    else:
        return path[0]   
	
#print("Waiting 1 minute before running")
#time.sleep(60)
data_dir = config["data_dir"]
if config["move"] is True:
    sample_sheet_path = config["samples"]
    sample_table = pd.read_csv(config["samples"], sep="\t", converters={'Experiment': str, "SampleID": str}).set_index("LabID")
    sample_table["path"] = sample_table.apply(lambda x: find_path(x.Experiment, x.SampleID, x.Barcode, data_dir), axis = 1)
    SAMPLES = sample_table.index.values

else:
    sample_table = pd.read_csv(config["samples"], sep="\t").set_index("LabID")
    sample_sheet_path = config["samples"]
    sample_table["path"] = sample_table.apply(lambda x: f"{data_dir}/{x.Experiment}/{x.SampleID}.fastq", axis=1)
    SAMPLES = sample_table.index.values  

SAMPLE_SHEET_PATH = sample_sheet_path
TIME = config["time"]# move to config file
SLEEP = "false"

print(SAMPLES, TIME)

include: "/workflow/rules/move_files.smk"
include: "/workflow/rules/centrifuge.smk"
include: "/workflow/rules/host_remove.smk"
include: "/workflow/rules/amr.smk"
include: "/workflow/rules/qc.smk"
include: "/workflow/rules/report.smk"
include: "/workflow/rules/cleanup.smk"
include: "/workflow/rules/mlst.smk"


rule all:
    input:
        expand("/mnt/reports/{sample}/{sample}_{time}_hours_report.pdf", sample = SAMPLES, time = TIME, sleep = SLEEP, sample_sheet_path = SAMPLE_SHEET_PATH),
        expand("/mnt/results/{sample}/{time}_hours/files/sample_sheet.tsv", sample = SAMPLES, time = TIME, sleep = SLEEP, sample_sheet_path = SAMPLE_SHEET_PATH),
        expand("/mnt/reports/{sample}/.{sample}_{time}_hours_report.pdf.tmp", sample=SAMPLES, time=TIME, sleep = SLEEP),
        expand("/mnt/results/{sample}/{time}_hours/microbial/{sample}_{time}_hours_hg38_removed.fastq.gz", sample=SAMPLES, time=TIME, sleep = SLEEP),
        expand("/mnt/results/{sample}/{time}_hours/centrifuge/centrifuge_raw.tsv.gz", sample=SAMPLES, time=TIME, sleep = SLEEP),
        #expand("results/{sample}/{time}_hours/transfer/transferred.txt", sample = SAMPLES, time = TIME),
        #expand("results/{sample}/{time}_hours/viral/viral_target_report.tsv", sample = SAMPLES, time = TIME),
        #expand("results/{sample}/{time}_hours/mlst/", sample= SAMPLES, time= TIME),
        #expand("results/{sample}/{time}_hours/amr/virulence_factor_summary.tsv", sample=SAMPLES, time=TIME)
