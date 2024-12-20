import pandas as pd


# SNAKEMAKE ARGUMENTS
# NODES = snakemake.config["taxonomy"]["nodes"]
# NAMES = snakemake.config["taxonomy"]["names"]
VIRULENCE_FILE = snakemake.input.vdf
SUMMARY_OUTPUT = snakemake.output.res


def gene_summary(df):
    gene_counts = df["GENE"].value_counts().to_dict()
    res = pd.DataFrame.from_dict(gene_counts, orient="index").reset_index()
    res.columns = ["Gene", "Counts"]
    return res

def main():
    df = pd.read_csv(VIRULENCE_FILE, sep="\t")
    if df.shape[0] != 0:
        res = gene_summary(df)
    else:
        # Failed sample
        res = pd.DataFrame({"Gene": ["N/A"], "Counts": ["N/A"]})
    res.to_csv(SUMMARY_OUTPUT, sep="\t", index=None)
    

if __name__ == "__main__":
    main()