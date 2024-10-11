import pandas as pd
import mappy as mp
import pyfastx
import os
import time
import taxonomy
from collections import Counter
import json
import traceback
from prettytable import PrettyTable

# SNAKEMAKE ARGUMENTS
NODES = snakemake.config["taxonomy"]["nodes"]
NAMES = snakemake.config["taxonomy"]["names"]
CENTRIFUGE_FILE = snakemake.input.file
CFG_THRESHOLD = snakemake.config["cfg_score"]
#VIRAL_THRESHOLD = snakemake.config["viral_score"]
TAXMETA = snakemake.config["taxonomy"]["speciesTaxMeta"]
FILENAMES = snakemake.config["taxonomy"]["speciesFileNames"]
FASTQ = snakemake.input.fastq
GENOME_DIR = snakemake.config["taxonomy"]["refseqDir"]
MULTI_OUTPUT = snakemake.output.multi
FAILED_OUTPUT = snakemake.output.failed
READ_OUTPUT = snakemake.output.read
REPORT_OUTPUT = snakemake.output.report
DICT_FILE = snakemake.config["taxonomy"]["dictFile"]
#VIRAL_TAXID = snakemake.config["taxonomy"]["viral_taxid"]
#SKIP_TAXID = snakemake.config["taxonomy"]["skip_taxid"]
THRESHOLD_CONFIG = snakemake.config["taxonomy"]["threshold_config"]
REPLACEMENT_LIST = snakemake.config["taxonomy"]["replacement_list"]


def remove_low_qual_reads(df, default_threshold, taxid_thresholds_path):
    """
    Filters reads with scores below the designated threshold, applying specific thresholds for taxids
    listed in the JSON files specified in the ranking.txt file. A threshold of '0' effectively skips
    the threshold check for specific taxids.

    Parameters
    ----------
    df : pandas.DataFrame
        Centrifuge output as pandas DataFrame.
    default_threshold : int
        Default centrifuge score threshold for false positive reads.
    taxid_thresholds_path : str
        Path to the directory containing the ranking.txt and JSON files.

    Returns
    -------
    tuple
        A tuple containing two elements:
        - above_threshold (pandas.DataFrame): Centrifuge dataframe with reads above score threshold.
        - failed (dict): Dictionary of failed reads due to low score.
          Format: {readID: {alignment: False, "reason": "score below threshold"}}.
    """
    
    print(f"Default threshold set to: {default_threshold}")

    # Read the ranking.txt file to determine the order of JSON files
    with open(os.path.join(taxid_thresholds_path, 'ranking.txt'), 'r') as f:
        json_files = [line.strip() for line in f]

    # Initialize an empty dictionary for taxid to score mapping
    taxid_score_dict = {}
    pretty_table = PrettyTable()
    pretty_table.field_names = ["Name", "Score"]

    # Read each JSON file as per the hierarchy and populate the dictionary
    for json_file in json_files:
        with open(os.path.join(taxid_thresholds_path, json_file), 'r') as f:
            data = json.load(f)
            name = data['Name']
            score = data['Score']
            taxids = data['TaxIDs']

            # Update the dictionary with taxid:score pairs
            for taxid in taxids:
                taxid_score_dict[taxid] = score

            # Add entry to the pretty table
            pretty_table.add_row([name, score])

    # Print the pretty table
    print(pretty_table)

    # Determine the appropriate threshold for each read
    def determine_threshold(taxid):
        return taxid_score_dict.get(taxid, default_threshold)

    # Apply filtering based on thresholds
    def apply_filtering(row):
        threshold = determine_threshold(row['taxID'])
        if threshold == 0:
            return True
        return row['score'] > threshold

    df['pass_filter'] = df.apply(apply_filtering, axis=1)

    # Separate above and below threshold reads
    above_threshold = df[df['pass_filter']].copy()
    below_threshold_df = df[~df['pass_filter']].copy()

    # Construct failed dictionary
    rejected_reads = below_threshold_df["readID"].unique()
    failed = {read_id: {"alignment": False, "reason": "score below threshold"} for read_id in rejected_reads}

    above_threshold = above_threshold.drop(columns=['pass_filter'])
    below_threshold_df = below_threshold_df.drop(columns=['pass_filter'])

    return above_threshold, failed

def species_id(taxID, tax, rank = "species"):
    """
    checks if taxID is a below taxonomy rank specified and returns that rank.

    Parameters
    ----------
    taxID
    tax
    rank

    Returns
    -------

    """
    taxobj = tax.parent(str(taxID), at_rank = rank)
    if taxobj is None:
        return tax.node(str(taxID))
    if taxobj.rank == rank:
        return taxobj
    return None


def get_name(taxID, tax):
    name = tax.node(str(taxID))
    if name is None:
        name = "NA"
        return name
    else:
        return name.name


def split_hits(df, tax):
    """
    Split centrifuge output into two categories: Unique hits and multi-match hits.
    Multi-match hits are further checked to see if hits are to the same taxonomy id/rank
    (ranks below species bought up to species)
    Also filters out reads below threshold
    Parameters
    ----------
    df
        centrifuge output as pandas DataFrame object
    tax
        taxonomy object from module taxonomy
    Returns
    -------
    unique_dict
        dictionary in format;  {readID: taxID}
    multi_dict
        dictionary in format; {readID:[taxIDs]}
    """

    # Split hits
    unique_df = df[df["numMatches"] == 1]
    multi_df = df[df["numMatches"] != 1]

    unique_dict = dict()
    multi_dict = dict()

    zipped = zip(unique_df["readID"], unique_df["taxID"])
    for r, t in zipped:
        unique_dict[r] = species_id(str(t), tax)

    # Check if multiple hits for read belong to same taxID i.e. Are actually unique hit
    multi_group = multi_df.groupby("readID")
    for read, dat in multi_group:
        tax_nodes = [species_id(str(t),tax) for t in dat["taxID"].unique() if species_id(str(t),tax) is not None]
        tax_ids = set(t.id for t in tax_nodes)
        if len(tax_ids) == 1:
            unique_dict[read] = tax_nodes[0]
        else:
            multi_dict[read] = tax_nodes

    # Bring taxIDS with ranks below species up to level

    return unique_dict, multi_dict


def create_fastq_dict(fastq):
    fastq_dict = {}
    for name, seq, qual, comment in pyfastx.Fastx(fastq):
        fastq_dict[name] = seq
    return fastq_dict


def create_tax_dict(dict_file):
    tax_dict = dict()
    df = pd.read_csv(dict_file, sep="\t")
    for index, row in df.iterrows():
        tax_dict[row["species_taxid"]] = row["filename"]
    return tax_dict


def make_aligner(tax_dict, genome_dir, tax_id):
    if tax_id not in tax_dict.keys():
        return None
    path = os.path.join(genome_dir, tax_dict[tax_id])
    if os.path.isfile(path):
        aligner = mp.Aligner(path, preset="map_ont", best_n=1)
        return aligner
    else:
        return None


def get_fastq(read_id, fastq_dict):
    fastq = fastq_dict[read_id]
    return fastq


def alignment(fastq, aligner):
    mapping = aligner.map(fastq)
    try:
        hit = next(mapping)
        res = {"hit": hit.ctg,
               "alignment_length": hit.blen,
               "total_matches": hit.mlen,
               "mapping_quality": hit.mapq,
               "identity": round(hit.mlen / hit.blen, 3),
               "alignment": True}
    except StopIteration:
        res = {"alignment": False,
               "reason": "No alignment found"}
    return res


def get_read_alignments(tax_ids, genome_dir, tax_dict, fastq, tax):
    hit_dict = dict()
    for tax_id in tax_ids:
        node = tax.node(str(tax_id))
        if node.rank == "species":
            tax_id = int(tax_id)
            aligner = make_aligner(tax_dict, genome_dir, tax_id)
            if aligner is None:
                hit_dict[tax_id] = {"alignment": False,
                                    "reason": "No reference genome"}
            else:
                hit_dict[tax_id] = alignment(fastq, aligner)
        else:
            hit_dict[tax_id] = {"alignment": False,
                                "reason": "tax above species level: {}".format(node.rank)}
    return hit_dict


def pick_best_alignment(align_dict):
    top = None
    best_identity = 0
    for k, v in align_dict.items():
        if v["alignment"] is True:
            if v["identity"] > best_identity:
                top = k
                best_identity = v["identity"]
    if top is None:
        return False
    else:
        return top

try:
    def main():
        start = time.time()
        print("loading data")
        tax = taxonomy.Taxonomy.from_ncbi(NODES, NAMES)
        df = pd.read_csv(CENTRIFUGE_FILE, sep="\t")
        # taxDict = create_tax_dict(TAXMETA, FILENAMES)
        taxDict = create_tax_dict(DICT_FILE)
        fastq_dict = create_fastq_dict(FASTQ)

        print("data loaded.... starting counts")

        above_threshold, failed = remove_low_qual_reads(df, CFG_THRESHOLD, THRESHOLD_CONFIG)
        report_dict, multi_hit_dict = split_hits(above_threshold, tax)

        print("Starting multi-match resolution")
        # Multi match sort
        multi_res = dict()
        for key, value in multi_hit_dict.items():
            fastq = get_fastq(key, fastq_dict)
            taxIDs = [int(t.id) for t in value]
            try:
                multi_res[key] = get_read_alignments(taxIDs, GENOME_DIR, taxDict, fastq, tax)
            except KeyError:
                failed[key] = taxIDs

        for k, v in multi_res.items():
            res = pick_best_alignment(v)
            if res:
                report_dict[k] = tax.node(str(res))
            else:
                failed[k] = v

        print("writing files")

        # Mutli match read assignments
        with open(MULTI_OUTPUT, "w") as f:
            json.dump(multi_res, f)

        # Failed reads
        with open(FAILED_OUTPUT, "w") as f:
            json.dump(failed, f)

        # Read assignments
        read_df = pd.DataFrame.from_dict(report_dict, orient="index").reset_index()
        read_df.columns = ["readID", "TaxID"]
        read_df["TaxID"] = read_df["TaxID"].apply(lambda x: x.id if (x is not None) else "NA")
        read_df["Organism"] = read_df["TaxID"].apply(lambda x: get_name(x,tax)) ## swap around
        read_df.to_csv(READ_OUTPUT, index=False, sep="\t")

        report_values = Counter([int(i.id) if i is not None else "No ID" for i in report_dict.values()])
        del report_values[9606]
        total_counts = sum(report_values.values())

        #Report output
        report_df = pd.DataFrame.from_dict(report_values, orient="index").reset_index()
        report_df.columns = ["Tax_ID", "Counts"]
        report_df["Organism"] = report_df["Tax_ID"].apply(lambda x: tax.node(str(x)).name if (tax.node(str(x)) is not None) else "ARGOS ISOLATE")
        report_df["Percentage"] = round(report_df["Counts"] / total_counts * 100, 3)
        report_df = report_df.sort_values(by="Percentage", ascending=False)
        report_df = report_df[["Organism", "Tax_ID", "Counts", "Percentage"]]

        # Path to CSV file in ref directory containing specific species/strain name replacements.
        species_complex_csv = pd.read_csv(REPLACEMENT_LIST)

        # Replace directly using the CSV data
        report_df['Organism'] = report_df['Organism'].replace(species_complex_csv.set_index('Original')['Replacement'])
        report_df.to_csv(REPORT_OUTPUT, index=False, sep="\t")


        end = time.time()
        print("Time to run is {} seconds".format(end - start))
except Exception as e:
    log_directory = '/mnt/'  # Update this path
    log_path = f'{log_directory}/error.log'
    with open(log_path, 'w') as f:
        f.write(traceback.format_exc())
        
if __name__ == "__main__":
    main()