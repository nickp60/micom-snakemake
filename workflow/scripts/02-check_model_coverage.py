import os
import sys
import argparse
import micom.qiime_formats as qf
from micom.qiime_formats import load_qiime_medium
from micom.workflows import grow, save_results, build, tradeoff, load_results
import seaborn as sns
import pandas as pd

def getargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('--abundances', required=True)
    parser.add_argument('--agora', required=True)
    parser.add_argument('--overrides')
    parser.add_argument('--taxrank', choices=["species", "genus"], default="species")
    parser.add_argument('--outdir', required=True)
    return  parser.parse_args()


if __name__ == "__main__":
    args = getargs()
    os.makedirs(args.outdir)
    taxrank = args.taxrank
    print("loading abundances")
    data = pd.read_csv(args.abundances)
    for col in ["sample_id",  taxrank, "abundance"]:
        assert col in data.columns, f"{col} must be a column name; colnames are" + " ,".join(data.columns)
    data = data[data.abundance > 0]
    # drop extra cols that could break the groupby summarizing later
    data = data[["sample_id", taxrank, "abundance"]]


    print("loading models")
    agora_taxa = qf.load_qiime_manifest(args.agora)

    print(agora_taxa.shape)
    # genera or species
    cols_of_interest = [taxrank]
    if taxrank == "species":
        cols_of_interest.append("id")
    agora_genera = agora_taxa[cols_of_interest]

    print("finding taxa in abundance file not present in models")
    ambig_genera_df = data[~data[taxrank].isin(agora_genera[taxrank])].query('abundance > 0').groupby(taxrank)["abundance"].agg(["median", "mean", "max", "count"]).sort_values("count", ascending=False)

    ambig_genera_df.to_csv(os.path.join(args.outdir, "unmatchable_taxa.csv"))
    agora_taxa[["id", taxrank]].to_csv(os.path.join(args.outdir, "available_taxa.csv"))

    if args.overrides:
        print("merging manual taxonomic overrides")
        replacementsdf = pd.read_csv(args.overrides, usecols=['A', 'B'])
        replacements = dict(zip(replacementsdf.A, replacementsdf.B))
        print(replacements)
        for k,v in replacements.items():
            tochange = data[taxrank] == k
            data.loc[tochange, taxrank] = v
    # especially after replacing taxa (but either way), we need to resum by taxrank to avoid duplicate rows
    print(f"recalculating abundance at {taxrank}")
    data = data.groupby(['sample_id', taxrank])['abundance'].sum().reset_index()


    print("checking for samples with no evaluable taxa")
    for sample in data.sample_id.unique():
        _data = data[data.sample_id == sample].query("abundance > 0")
        ntax = _data.shape[0]
        _data = _data[_data[taxrank].isin(agora_genera[taxrank])]

        if _data.shape[0] == 0:
            print(f"sample {sample} has no evaluable genera")
            data = data[data.sample_id != sample]
        else:
            print(f"sample {sample} has {_data.shape[0]}/{ntax} evaluable genera")

    # remove non-agora taxa
    print("generating found fraction plot")
    data = data[data[taxrank].isin(agora_genera[taxrank])]
    data_summary = data.groupby('sample_id').agg({'abundance': 'sum'})
    data_summary["tmp"] = "1"
    p = sns.boxplot(data=data_summary, x="abundance", y="tmp")
    sns.swarmplot(data=data_summary, x="abundance",  y="tmp", ax=p)
    p.figure.savefig(os.path.join(args.outdir, "MICOM_found_abundance_fraction.png"))
    print("saving workflow input table")
    data.to_csv(os.path.join(args.outdir, "postQC_abundances.csv"))
