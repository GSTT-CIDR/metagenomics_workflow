import pandas as pd
import taxonomy

# SNAKEMAKE ARGUMENTS
NODES = snakemake.config["taxonomy"]["nodes"]
NAMES = snakemake.config["taxonomy"]["names"]
REPORT = snakemake.input.report
BAC_REPORT = snakemake.output.bacterial
VIRAL_REPORT = snakemake.output.viral

def get_superkingdom(taxID, tax, rank = "superkingdom"):
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
        return None
    if taxobj.rank == rank:
        return str(taxobj.name)
    return None

def main():
    tax = taxonomy.Taxonomy.from_ncbi(NODES, NAMES)
    df = pd.read_csv(REPORT, sep="\t")
    df["superkingdom"] = df["Tax_ID"].apply(lambda x: get_superkingdom(x, tax))

    # Create a copy when slicing the DataFrame to avoid setting with copy warning
    viral = df[df["superkingdom"] == "Viruses"].copy()
    viral_counts = viral["Counts"].sum()
    viral["Percentage"] = viral["Counts"] / viral_counts * 100

    bac = df[df["superkingdom"] != "Viruses"].copy()
    bac_counts = bac["Counts"].sum()
    bac["Percentage"] = bac["Counts"] / bac_counts * 100

    viral.to_csv(VIRAL_REPORT, sep="\t", index=None)
    bac.to_csv(BAC_REPORT, sep="\t", index=None)


if __name__ == "__main__":
    main()
