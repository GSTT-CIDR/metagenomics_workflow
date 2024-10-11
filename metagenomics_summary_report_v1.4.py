import pandas as pd
import glob
import os
from collections import Counter
import re
import taxonomy
import csv
import argparse
from styleframe import StyleFrame

# Setup argparse
parser = argparse.ArgumentParser(description='The script produces a summary of outputs from the metagenomics pipeline in XLSX or format. Sample names should be provided to the script one per line and must be identical to the sample names in the sample_sheet used to launch the metagnomics pipeline.')
parser.add_argument('--sample_list', type=str, help='File path to the list of sample names.')
parser.add_argument('--results_directory', type=str, help='Directory path for the results.')
parser.add_argument('--output_directory', type=str, help='Output file path.')
parser.add_argument('--abundance_threshold', type=float, default=1.0, help='Minimum abundance threshold for filtering (default: 1.0).')
parser.add_argument('--delimiter', type=str, default='False', help='Set the list delimiter to \';\' - helps when manipulating outside of Excel.')
parser.add_argument('--sort_by_date', type=str, default='False', help='Detect YY/MM/DD date within sample name and sort by it in ascending order.')
parser.add_argument('--as_xlsx', type=str, default='False', help='Export as XLSX instead of CSV')
args = parser.parse_args()


# Replace hardcoded paths with arguments
NODES = "/mnt/db/ref/refseq/taxonomy/nodes.dmp"
NAMES = "/mnt/db/ref/refseq/taxonomy/names.dmp"
tax = taxonomy.Taxonomy.from_ncbi(NODES, NAMES)
### TIMEPOINTS TO EXTRACT DATA FROM
time_points = [0.5, 1, 2, 16, 24]

def main():
    # Reading sample names from the file specified in the command line
    sample_names = read_sample_names(args.sample_list)
    print("Samples specified for analysis:")
    print(sample_names)

    # Use results_directory from arguments
    sample_paths = [os.path.join(args.results_directory, sample) for sample in sample_names]
    existing, non_existing = check_directories(sample_paths)

    files = existing
    print("Samples not present in target directory:")
    print(non_existing)

    # Use abundance_threshold from arguments
    threshold = args.abundance_threshold

    results = create_df(files, time_points, threshold)

    # Reset the index to make the sample names a regular column instead of the index
    master = pd.DataFrame.from_dict(results, orient="index").reset_index()
    
    # Rename the new column to "Sample"
    master.rename(columns={'index': 'Sample'}, inplace=True)    
    
    # Save the dataframe in the specified format
    export(master, args)


def read_sample_names(sample_list_file):
    sample_names = []
    with open(sample_list_file, 'r') as file:
        for line in file:
            sample_names.append(line.strip())
    return sample_names

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    return [atoi(c) for c in re.split(r'(\d+)', text)]

def check_directories(dir_list):
    existing_dirs = []
    non_existing_dirs = []
    for directory in dir_list:
        if os.path.isdir(directory):
            existing_dirs.append(directory)
        else:
            non_existing_dirs.append(directory)
    return existing_dirs, non_existing_dirs

# Function to parse the NanoStat summary file
def parse_nanostat_summary(filename):
    metrics = {}
    with open(filename, 'r') as file:
        content = file.read()
        # Extracting specific metrics
        q15_matches = re.search(r'>Q15:\s+\d+ \(([\d\.]+)%\)', content)
        median_quality_matches = re.search(r'Median read quality:\s+([\d\.]+)', content)
        total_bases_matches = re.search(r'Total bases:\s+([\d\.,]+)', content)
        n50_matches = re.search(r'Read length N50:\s+([\d\.,]+)', content)
        
        if q15_matches:
            metrics['Q15_Reads'] = float(q15_matches.group(1))
        if median_quality_matches:
            metrics['Median_Read_Quality'] = float(median_quality_matches.group(1))
        if total_bases_matches:
            metrics['Total_Bases'] = float(total_bases_matches.group(1).replace(',', ''))
        if n50_matches:
            metrics['Sequencing_N50'] = float(n50_matches.group(1).replace(',', ''))
    
    return metrics

def create_df(files, time_points, threshold):
    results_dir = args.results_directory
    results = dict()
    for f in files:
        sample = f.split("/")[-1]
        sample_sheet_path = f"{results_dir}/{sample}/0.5_hours/files/sample_sheet.tsv"
        
        sample_info_dict = {}  # Initialize an empty dictionary
        
        if os.path.exists(sample_sheet_path):
            sample_sheet_df = pd.read_csv(sample_sheet_path, sep="\t")
            sample_info = sample_sheet_df[sample_sheet_df['LabID'] == sample]
            
            if not sample_info.empty:
                sample_info_dict = sample_info.to_dict(orient='records')[0]  # Convert the row to a dictionary
            else:
                sample_info_dict = {
                    'Experiment': 'sample sheet not found',
                    'SampleID': 'sample sheet not found',
                    'Barcode': 'sample sheet not found',
                    'AnonymisedIdentifier': 'sample sheet not found',
                    'CollectionDate': 'sample sheet not found',
                    'SampleClass': 'sample sheet not found',
                    'SampleType': 'sample sheet not found',
                    'Operator': 'sample sheet not found',
                    'Notes': 'sample sheet not found'
                }
        else:
            sample_info_dict = {
                'Experiment': 'sample sheet not found',
                'SampleID': 'sample sheet not found',
                'Barcode': 'sample sheet not found',
                'AnonymisedIdentifier': 'sample sheet not found',
                'CollectionDate': 'sample sheet not found',
                'SampleClass': 'sample sheet not found',
                'SampleType': 'sample sheet not found',
                'Operator': 'sample sheet not found',
                'Notes': 'sample sheet not found'
            }

        # Safely retrieve fields with .get(), providing a default value if the field is missing
        sample_info_dict = {
            'Experiment': sample_info_dict.get('Experiment', 'not available'),
            'SampleID': sample_info_dict.get('SampleID', 'not available'),
            'Barcode': sample_info_dict.get('Barcode', 'not available'),
            'AnonymisedIdentifier': sample_info_dict.get('AnonymisedIdentifier', 'not available'),
            'CollectionDate': sample_info_dict.get('CollectionDate', 'not available'),
            'SampleClass': sample_info_dict.get('SampleClass', 'not available'),
            'SampleType': sample_info_dict.get('SampleType', 'not available'),
            'Operator': sample_info_dict.get('Operator', 'not available'),
            'Notes': sample_info_dict.get('Notes', 'not available')
        }

        results[sample] = sample_info_dict

        for time in time_points:
            bac_path = f"{f}/{time}_hours/centrifuge/bacterial_centrifuge_report.tsv"
            read_stats_path = f"{f}/{time}_hours/host/{sample}_{time}_hours_map_stats.txt"
            qc_stats_path = f"{f}/{time}_hours/qc/nanostat_summary.txt"
            
            if os.path.exists(read_stats_path):
                with open(read_stats_path, 'r') as stats_file:
                    total_reads = int(stats_file.readline().strip())
                    human_reads = int(stats_file.readline().split('/')[0].strip())
            else:
                total_reads, human_reads = 0, 0
            
            if total_reads > 0:
                human_reads_percentage = round((human_reads / total_reads) * 100, 1)
            else:
                human_reads_percentage = 'file not present'
            
            if os.path.exists(qc_stats_path):
                qc_metrics = parse_nanostat_summary(qc_stats_path)
            else:
                qc_metrics = {'Q15_Reads': 0, 'Median_Read_Quality': 0.0, 'Total_Bases': 0, 'Sequencing_N50': 0}

            if os.path.exists(bac_path):
                df = pd.read_csv(bac_path, sep="\t")
                if df.empty:
                    total_counts, organisms, counts, percentage = 0, 0, 0, 0
                else:
                    total_counts = df["Counts"].sum()
                    df["Percentage"] = round(df["Counts"] / total_counts * 100, 3)
                    df = df[(df["Percentage"] >= threshold) | (df["Organism"].str.contains("Candida|Aspergillus|Mycoplasma|Legionella|Chlamydia|Pneumocystis|Mycobacterium|Mycobacteroides"))]
                    df = df[["Organism", "Counts", "Percentage"]]
                    organisms, counts, percentage = df.apply(lambda x: "\n".join([str(i) for i in x]))
            else:
                total_counts, organisms, counts, percentage = 'file not present', 'file not present', 'file not present', 'file not present'
            
            read_path = f"{f}/{time}_hours/centrifuge/viral_centrifuge_report.tsv"
            virus, v_counts = 'file not present', 'file not present'
            
            if os.path.exists(read_path):
                read_df = pd.read_csv(read_path, sep="\t")
                if read_df.empty:
                    virus, v_counts = 0, 0
                else:
                    read_df = read_df[["Organism", "Counts"]]
                    virus, v_counts = read_df.apply(lambda x: "\n".join([str(i) for i in x]))
            
            results[sample].update({
                f"Total reads {time} hrs": total_reads,
                f"Human reads {time} hrs": human_reads,
                f"Human reads (%) {time} hrs": human_reads_percentage,
                f"Total classified reads {time} hrs": total_counts,
                f"Sequencing N50 (bp) {time} hrs": qc_metrics['Sequencing_N50'],
                f"Proportion >Q15 quality (%) {time} hrs": qc_metrics['Q15_Reads'],
                f"Median read quality (PHRED score) {time} hrs": qc_metrics['Median_Read_Quality'],
                f"Total bases (bp) {time} hrs": qc_metrics['Total_Bases'],
                f"Organisms (excluding viruses) {time} hrs": organisms,
                f"Organisms (excluding viruses) read counts {time} hrs": counts,
                f"Organism (excluding viruses) percentage abundance {time} hrs": percentage,
                f"Viral organisms {time} hrs": virus,
                f"Viral read counts {time} hrs": v_counts
            })
    return results

def replace_delimiter_if_true(df, args):
    """
    Replace '\n' with ';' in the DataFrame if args.delimiter is 'True'.
    
    Parameters:
    df (pd.DataFrame): The DataFrame containing cells with delimited lists.
    args: An object with a 'delimiter' attribute.
    
    Returns:
    pd.DataFrame: The DataFrame with '\n' replaced by ';' if args.delimiter is 'True'.
    """
    if args.delimiter == 'True':
        df = df.applymap(lambda x: x.replace('\n', ';') if isinstance(x, str) else x)
    return df

def export(master, args):
    if args.as_xlsx == 'True':
        # Replacing line for StyleFrame 
        master_xlsx = master.replace('\r', '\n')
        # Write to Excel with StyleFrame for correct line breaks
        StyleFrame(master_xlsx).to_excel(args.output_directory).save()
    else: 
        # Export to CSV - with added delim function
        master_export = replace_delimiter_if_true(master, args)
        master_export.to_csv(args.output_directory, index=False)



# Entry point of the script
if __name__ == "__main__":
    main()
