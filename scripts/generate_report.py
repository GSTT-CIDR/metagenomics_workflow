import pandas as pd
import csv
from jinja2 import Environment, FileSystemLoader
import re
from weasyprint import HTML, CSS
from datetime import datetime

SAMPLE = snakemake.wildcards.sample
INTERVAL = float(snakemake.wildcards.time)
CFG_PATH = snakemake.input.centrifuge
CFG_RAW_PATH = snakemake.input.centrifuge_raw
VIRAL_PATH = snakemake.input.viral
AMR_SUMMARY = snakemake.input.amr_summary
AMR_REPORT = snakemake.input.amr_report
VF_PATH = snakemake.input.vf
QC_PATH = snakemake.input.qc
OUTPUT = str(snakemake.output)
REPORT_HTML = snakemake.config["pdf"]["html"]
REPORT_CSS = snakemake.config["pdf"]["css"]
BOOTSTRAP_CSS = snakemake.config["pdf"]["bootstrap"]
SAMPLE_TABLE = snakemake.config["samples"]
SAMTOOLS_STAT = snakemake.input.stats
THRESHOLD = snakemake.config["abundance_threshold"]
CFG_THRESHOLD = snakemake.config["cfg_score"]
TARGETS = snakemake.config["targets"]
VERSION = snakemake.config["pdf"]["version"]



def convert_bp(size):
    size = float(str(size).replace(",", ""))
    for x in ["bp", "Kb", "Mb", "Gb", "Tb"]:
        if size < 1000.0:
            if x == "bp":
                return "{:.0f} {}".format(size, x)
            else:
                return "{:.2f} {}".format(size, x)
        size /= 1000.0
    return size


def summary_qc(path):
    qc_dict = dict()
    with open(path, "r") as f:
        data = f.read().splitlines()
    for line in data:
        if line.startswith("Mean read length"):
            dat = re.split(r"\s{2,}", line)[1]
            qc_dict["mean_read_len"] = dat
        elif line.startswith("Mean read quality"):
            dat = re.split(r"\s{2,}", line)[1]
            qc_dict["mean_read_qual"] = dat
        elif line.startswith("Number of reads"):
            dat = re.split(r"\s{2,}", line)[1]
            qc_dict["nano_reads"] = dat
        elif line.startswith("Total bases"):
            dat = re.split(r"\s{2,}", line)[1]
            qc_dict["total_bp"] = convert_bp(dat)
    return qc_dict


def samtools_stats(path):
    with open(path, "r") as f:
        data = f.read().splitlines()
    total_reads = data[0].strip()  # Reads the first line for total read count
    human_reads = data[1].strip()  # Reads the second line for human read count
    stats = {"total_reads": total_reads, "human_reads": human_reads}
    return stats



def patient_info(path, id):
    df = pd.read_csv(path, sep = "\t")
    df = df[df["LabID"] == str(id)]
    sample_dict = {"LabID" : df["LabID"].values[0],
                   "Experiment": df["Experiment"].values[0],
                   "SampleID": df["SampleID"].values[0],
                   "Barcode": df["Barcode"].values[0],
                   "Version": VERSION,
                   "SampleType": df["SampleType"].values[0],
                   "SampleClass": df["SampleClass"].values[0],
                   "Operator": df["Operator"].values[0],
                   "Notes": df["Notes"].values[0],
                   }
    return sample_dict

# def is_target(s, target_file = TARGETS):
#     with open(TARGETS, "r") as f:
#         targets = [t.strip() for t in f]
#     if s["Organism"] in targets:
#         return ["color: red"] * 3
#     else:
#         return ["color: black"] * 3

# Change this from hardcoding threshold
def cfg_to_html(path, threshold = THRESHOLD):
    exceptions = "Aspergillus|Candida|Chlamydia|Pneumocystis|Mycoplasma"
    cfg_dict = dict()
    df = pd.read_csv(path, sep="\t")
    df = df.round(2)
    df = df.drop(columns=["Tax_ID", "superkingdom"])
    cfg_dict["cfg_full"] = df.to_html(classes="table table-striped", border=0, justify="left", index=False)
    cfg_dict["micro_reads"] = sum(df["Counts"])
    ic = df[df["Organism"] == "Jonesia denitrificans"]
    if ic.empty:
        cfg_dict["ic"] = "NA/NA"
    else:
        cfg_dict["ic"] = "{}/{}%".format(int(ic["Counts"]), float(ic["Percentage"]))
    above_df = df.copy()
    above_df = above_df[(above_df["Percentage"] > threshold) |(above_df["Organism"].str.contains(exceptions))]
    cfg_dict["cfg_top"] = above_df.to_html(classes="table table-striped", border=0, justify="left", index=False)
    return cfg_dict


def unclassified_reads(path, threshold=CFG_THRESHOLD):
    df = pd.read_csv(path, sep="\t")
    unclass = len(df[df["seqID"].str.contains("unclassified")]["readID"].unique())
    below = len(df[df["score"] < threshold]["readID"].unique()) - unclass
    return {"unclass_reads": unclass,
            "reads_below_threshold": below}



def amr_summary(path):
    summary_dict = dict()
    df = pd.read_csv(path, sep="\t", header=None, names=["Species", "Gene", "Hits"])
    res = df[df["Hits"] >= 2]
    if res.empty:
        res = res.append({"Species": "No genes above threshold", "Gene": "NA", "Hits": "NA"}, ignore_index=True)
    summary_dict["amr_summary"]  = res.to_html(classes="table table-striped", border=0, justify="left", index=False)
    return summary_dict

def amr_report(path, summary_path):
    report_dict = dict()
    with open(path, "r") as f:
        dat = list(filter(None, [i for i in csv.reader(f, delimiter="\t")]))
    column_names = ["SPECIES", "SEQUENCE", "START", "END", "GENE", "%COVERAGE"]
    master_df = pd.DataFrame(columns = column_names)
    if dat[0][0].startswith("#FILE"):
        with open(summary_path, "r") as f:
            species = [i for i in csv.reader(f, delimiter="\t")][0][0]
        df = pd.read_csv(path, sep="\t", usecols=["SEQUENCE", "START", "END", "GENE", "%COVERAGE"])
        df["SPECIES"] = species
        report_dict["amr_report"]  = df[column_names].to_html(classes="table table-striped", border=0, justify="left", index=False)
        return report_dict
    elif dat[0][0].startswith("Results for species"):
        coord = []
        species = []
        pos = []
        no_res = []
        for num, line in enumerate(dat):
            if line[0].startswith("Results for species"):
                species.append(line[1])
                if len(pos) == 1:
                    pos.append(num)
                    coord.append(pos)
                    pos = []
            if line[0].startswith("#FILE"):
                pos.append(num)
            if line[0].startswith("No results"):
                no_res.append(line[0])
        if len(pos) == 1:
            pos.append(len(dat))
            coord.append(pos)
        if len(coord) == 0:
            df = pd.DataFrame(columns=["SPECIES", "RESULTS"])
            for s, n in zip(species, no_res):
                df = df.append({"SPECIES": s,
                                "RESULTS": n}, ignore_index=True)
            report_dict["amr_report"]  = df.to_html(classes="table table-striped", border=0, justify="left", index=False)
            return report_dict
        for s, c in zip(species,coord):
            tmp = dat[c[0]:c[1]]
            df = pd.DataFrame(tmp[1:], columns=tmp[0])
            df["SPECIES"] = s
            master_df = master_df.append(df[column_names])
        report_dict["amr_report"]  = master_df.to_html(classes="table table-striped", border=0, justify="left", index=False)
    else:
        report_dict["amr_report"] = "No Results"
    return report_dict

def viral_report(path):
    viral_dict = dict()
    df = pd.read_csv(path, sep="\t", usecols=["Organism", "Counts"])
    ic = df[df["Organism"] == "Tobacco mosaic virus"]
    viral_dict["viral_reads"] = sum(df["Counts"])
    if ic.empty:
        viral_dict["ic"] = "NA/NA"
    else:
        viral_dict["ic"] = "{}".format(int(ic["Counts"]))
    viral_dict["viral_report"] = df.to_html(classes="table table-striped", border=0, justify="left", index=False)
    return viral_dict

def virulence_factors(path):
    vf_dict = dict()
    df = pd.read_csv(path, sep="\t")
    vf_dict["vf_report"] = df.to_html(classes="table table-striped", border=0, justify="left", index=False)
    return vf_dict

# Initialize report_dict with basic information
report_dict = {"time": str(INTERVAL) + " hrs",
               "title": "Clinical metagenomics report",
               "date": datetime.now()}

# Update the dictionary with various reports
report_dict.update(samtools_stats(SAMTOOLS_STAT))
report_dict.update(patient_info(SAMPLE_TABLE, SAMPLE))

# Get bacterial and viral reads, and update the dictionary
bacterial_data = cfg_to_html(CFG_PATH)
viral_data = viral_report(VIRAL_PATH)

report_dict.update(bacterial_data)
report_dict.update(viral_data)

# Calculate the sum of bacterial and viral reads and add it to the dictionary
# Ensure that 'micro_reads' and 'viral_reads' are integers before summing
sum_reads = bacterial_data.get('micro_reads', 0) + viral_data.get('viral_reads', 0)
report_dict['sum_reads'] = sum_reads

# Continue updating the dictionary with other reports
report_dict.update(summary_qc(QC_PATH))
report_dict.update(unclassified_reads(CFG_RAW_PATH))
report_dict.update(amr_summary(AMR_SUMMARY))
report_dict.update(amr_report(AMR_REPORT, AMR_SUMMARY))
report_dict.update(virulence_factors(VF_PATH))



env = Environment(loader=FileSystemLoader(".")) # Change to "." for grid
# declare our jinja template
template = env.get_template(REPORT_HTML)

html_out = template.render({"report": report_dict})


pdf_name = OUTPUT
HTML(string=html_out).write_pdf(pdf_name,
                                stylesheets=[CSS(REPORT_CSS),
                                             CSS(BOOTSTRAP_CSS)])
