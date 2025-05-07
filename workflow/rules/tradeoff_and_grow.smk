import os
import sys


from micom import Community
import micom.qiime_formats as qf
from micom.qiime_formats import load_qiime_medium
from micom.workflows import grow, save_results, build, tradeoff, load_results
from micom import Community
import seaborn as sns
import pandas as pd




data = pd.read_csv(config["abundances"])
for col in ["sample_id",  config['taxrank'], "abundance"]:
    assert col in data.columns, f"{col} must be a column name"
data = data[data.abundance > 0]
# drop extra cols that could break the groupby summarizing later
data = data[["sample_id", config['taxrank'], "abundance"]]

metadata = pd.DataFrame(data.sample_id.unique(), columns=["sample_id"])
if "ignore_samples" in config:
    if config["stage"] == "grow":
        print("Are you sure you want to exclude samples from the grow stage?")
    for sample in config["ignore_samples"]:
        print("dropping sample")
        metadata = metadata[metadata.sample_id != sample]

if "testn" in config:
    metadata = metadata.head(config["testn"])


agora_taxa = qf.load_qiime_manifest(config["agora"])
print(agora_taxa[config["taxrank"]])

#print(agora_taxa.shape)
#agora_taxa = qf.load_qiime_manifest("/lila/data/brinkvd/users/watersn/micom-databases/recipes/agora2/databases/agora201_gtdb220_species_1.qza")
#agora_taxa = qf.load_qiime_manifest("/lila/data/brinkvd/users/watersn/micom-databases/recipes/carveme/databases/carveme260d0f1_gtdb220_species_1.qza")
#agora_taxa = qf.load_qiime_manifest("../resources/agora201_gtdb207_species_1.qza")
print(agora_taxa.shape)
# genera or species
if config["taxrank"] == "species":
    cols_of_interest = ["id", config['taxrank']]
else:
    cols_of_interest = [config['taxrank']]
agora_genera = agora_taxa[cols_of_interest]

ambig_genera_df = data[~data[config['taxrank']].isin(agora_genera[config["taxrank"]])].query('abundance > 0').groupby(config["taxrank"])["abundance"].agg(["median", "mean", "max", "count"]).sort_values("count", ascending=False)

ambig_genera_df.to_csv("unmatchable_taxa.csv")
agora_taxa[["id", config['taxrank']]].to_csv("available_taxa.csv")
# def make_autokey(data, agora_genera, rank="species"):
#     uniq_taxa = pd.DataFrame(data[rank].unique(), columns=["raw"])
#     for i in in uniq_taxa.shape[0]:
#     	taxon = uniq_taxa.iloc[i, 0]
# 	genus  = taxon.split(" ")[0]
# 	species  = " ".join(taxon.split(" ")[1:])
# 	ref_genus =

if "taxakey" in config:
#    if config["taxakey"] == "auto":
#       data = make_autokey(data)
    replacementsdf = pd.read_csv(config["taxakey"], usecols=['A', 'B'])
    replacements = dict(zip(replacementsdf.A, replacementsdf.B))
    print(replacements)
    for k,v in replacements.items():
        tochange = data[config['taxrank']] == k
        data.loc[tochange, config['taxrank']] = v
# especially after replacing taxa (but either way), we need to resum by taxrank to avoid duplicate rows
data = data.groupby(['sample_id', config['taxrank']])['abundance'].sum().reset_index()


# check for samples with no evaluable taxa
for sample in data.sample_id.unique():
    _data = data[data.sample_id == sample].query("abundance > 0")
    ntax = _data.shape[0]
    _data = _data[_data[config['taxrank']].isin(agora_genera[config["taxrank"]])]

    if _data.shape[0] == 0:
        print(f"sample {sample} has no evaluable genera")
        data = data[data.sample_id != sample]
    else:
        print(f"sample {sample} has {_data.shape[0]}/{ntax} evaluable genera")

# remove non-agora taxa
data = data[data[config['taxrank']].isin(agora_genera[config["taxrank"]])]
data_summary = data.groupby('sample_id').agg({'abundance': 'sum'})
data_summary["tmp"] = "1"
p = sns.boxplot(data=data_summary, x="abundance", y="tmp")
sns.swarmplot(data=data_summary, x="abundance",  y="tmp", ax=p)
p.figure.savefig("MICOM_found_abundance_fraction.png")
assert config["stage"] in ["tradeoff", "grow"], "stage must either be tradeoff or grow"

if config["stage"] == "tradeoff":
    targets = [
        "tradeoff.csv",
        "tradeoff.html",
        ]
else:
    assert config["tradeoff"] is not None, "Must set tradeoff config value"
    targets = "growth.zip"

rule all:
    input:
        "manifest.csv",
        expand("models/{sample}.pickle", sample = metadata.sample_id),
        targets



rule build:
    input:
        agora_file=config["agora"],
    output:
        outf="models/{sample}.pickle",
        out_manifest = "manifests/{sample}.csv"
    threads: 2
    resources:
        mem_mb=lambda wc, attempt: 8 * 1024 * attempt,
    params:
        cutoff=config["cutoff"]
    run:
        from micom.logger import logger
        logger.setLevel("INFO")
        thisdata = data[data.sample_id == wildcards.sample]
        print(thisdata)
        manifest_osqp = build(
            thisdata,
            input.agora_file,
            os.path.dirname(output.outf),
            solver="osqp",
            cutoff=params.cutoff,
            threads=1)
        manifest_osqp.to_csv(output.out_manifest, index=False)



rule merge_minifests:
    input:
        infiles = expand("manifests/{sample}.csv", sample = metadata.sample_id)
    output:
        outf = "manifest.csv"

    run:
        manifest = pd.concat((pd.read_csv(f) for f in input.infiles), ignore_index=True)
        manifest.to_csv(output.outf, index=False)

rule get_tradeoff:
    input:
        manifest = "manifests/{sample}.csv",
        model = "models/{sample}.pickle",
        medium = config["medium"]
    output:
        outf = "tradeoffs/{sample}.csv"
    threads: 2
    resources:
        mem_mb=lambda wc, attempt: 8 * 1024 * attempt,
        runtime=lambda wc, attempt: 24 * 60 * attempt,
    run:
        medium = load_qiime_medium(input.medium)
        manifest = pd.read_csv(input.manifest)
        tradeoff_results = tradeoff(
            manifest,
            os.path.dirname(input.model),
            medium,
            threads=threads
        )
        tradeoff_results.to_csv(output.outf, index=False)

use rule merge_minifests as merge_tradeoffs with:
    input:
        infiles = expand("tradeoffs/{sample}.csv", sample = metadata.sample_id)
    output:
        outf = "tradeoff.csv"

rule plot_tradeoff:
    input:
        inf = "tradeoff.csv"
    output:
        html = "tradeoff.html"
    run:
        from micom.viz import plot_tradeoff
        tradeoff = pd.read_csv(input.inf)
        plot_tradeoff(tradeoff, filename=output.html,  tolerance=1e-06)

rule grow_sample:
    input:
        manifest = "manifests/{sample}.csv",
        model = "models/{sample}.pickle",
        medium = config["medium"]
    output:
        growth="growth/{sample}.zip"
    params:
        tradeoff = config["tradeoff"]
    threads: 4
    resources:
        runtime=lambda wc, attempt: 4 * 60 * attempt,
    run:
        medium = load_qiime_medium(input.medium)
        manifest = pd.read_csv(input.manifest)
        growth = grow(
            manifest,
            os.path.dirname(input.model),
            medium,
            tradeoff=params.tradeoff,
            threads=threads,
            presolve=True)
        save_results(growth, output.growth)

def merge_gr(gra, grb):
    from copy import deepcopy
    new = deepcopy(gra)
    _exchanges = pd.concat([gra.exchanges, grb.exchanges], ignore_index=True)
    _growth_rates =  pd.concat([gra.growth_rates, grb.growth_rates], ignore_index=True)
    _annotations =  pd.concat([gra.annotations, grb.annotations], ignore_index=True).drop_duplicates().reset_index(drop=True)
    new.exchanges = _exchanges
    new.growth_rates = _growth_rates
    new.annotations = _annotations
    return(new)


rule merge_growth:
    input:
        infiles = expand("growth/{sample}.zip", sample= metadata.sample_id)
    output:
        growth = "growth.zip"
    run:
        growth = load_results(input.infiles[0])
        for i in range(2, len(input.infiles)):
            print(i)
            growthb = load_results(input.infiles[i])
            growth = merge_gr(growth, growthb)
            print(growth.growth_rates.shape)
        save_results(growth, output.growth)
