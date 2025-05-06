# micom-snakemake
A snakemake workflow for running MICOM  https://micom-dev.github.io/micom/


It can be tricky to get your data formated for running MICOM, and several steps of MICOM can be quite time consuming. Luckily, they are also parallelizable! Here we describe how to use the included Snakemake workflow to distribute those tasks to the cloud or HPC, and hand the results back to the tutorial.



## prereqs

### A conda env:
```
mamba create --name micomwf -c conda-forge -c bioconda micom snakemake seaborn
```

Then add some deps easier to manage with pip:
```
#Biom-format requires CFLAGS="-std=c99"
export CFLAGS="-std=c99"
pip install numpy==2.1.1 Cython biom-format==2.1.16 micom seaborn
```

Be sure to also include an executor relevant to your environment; I'm using `snakemake-executor-plugin-cluster-generic`.

### Snakemake

### Singularity

## Prepare 16S amplicon data


After your sequences-by-counts table has been generated (DADA2, DEBLUR), extract the ASV/OTU sequences as a fasta file.  These can be annotated by blasting against a database created from the AGORA genomes.


### Download AGORA2 genomes
Download AGORA2 genomes as described [here](https://github.com/vdblab/resources/tree/main/workflow/sylphflux)

### extract 16S and build BLAST database

This step uses barrnap to pull each genome's rRNAs into a fasta file, adjusts the headers, merges, and creates a blast db under `<directory>/db/<dbname>`.

```
snakemake --config stage=format16S  agora_genomes=$PWD/.tests/genomes/ dbname="testdb"  --directory $PWD/mydb/  -n
```

### annotate your ASVs
This step takes in your counts, your asv fasta, and the manifest of the database you used (see https://raw.githubusercontent.com/micom-dev/databases/refs/heads/main/recipes/agora2/manifests/ for instance).  This allows us to match the genome IDs in the blast database with the model names in the micom model bundles:

```

snakemake --config stage=annotate16S  id=batch1 asv_counts=$PWD/.tests/asv_counts.csv asv_fasta=$PWD/.tests/testasvs.fasta  model_manifest=$PWD/agora2_refseq_species.tsv db_path_nsq=$PWD/mydb/db/testdb.nsq  --directory $PWD/results/
```

## Prepare shotgun metagenomic data

### Get a Shotgun Taxonomic Profile
We use sylph to get taxonomic profiles against a database constructed from the AGORA2 genomes. See the workflow [here](https://github.com/vdblab/resources/tree/main/workflow/sylphflux) for downloading the genomes and formatting into a sylph sketch.  Extract a table from the results that is in "tall" format, with columns for `sample_id`, `species/genus`, and `abundance`.


## Merge with names from model manifest
The genome names and the names used for the model are often dissimilar; this (horrible) script attempts to match up the names of the genome ids with the model IDs based on the model manifests found in the `micom-dev/databases/` repo.
```
wget https://raw.githubusercontent.com/micom-dev/databases/refs/heads/main/recipes/agora2/manifests/agora2_refseq_species.tsv
# can run this with system R if desired and packages are present
singularity exec -B $PWD  docker://rocker/tidyverse:4.1 Rscript $PWD/workflow/scripts/genomeid_to_modelid_amplicon.R $PWD/results/annotate/batch1_abundances.csv $PWD/agora2_refseq_species.tsv $PWD/results/annotate/batch1_abundances_renamed.csv
```


# Usage
Now that you have an abundances table either from 16S or shotgun data, its time to MICOM!

## Input Data

This workflow requires three files:

- a long-format sample abundances table with columns as follows:
  - `sample_id`
  - `genus`
  - `abundance`
- a MICOM-formatted AGORA model bundle: `wget -O agora201_refseq216_species_1.qza https://zenodo.org/records/7739096/files/agora201_refseq216_species_1.qza?download=11`
- a MICOM-formatted diet/medium: `wget -O western_diet_gut_agora.qza https://github.com/micom-dev/media/raw/main/media/western_diet_gut_agora.qza`


<!-- #- (optional) file containing genus taxonomy translations between the tool generating the counts and the underlying taxa models used by MICOM. Column `A` shoudl have the genus as provided by your tool, and column `B` should have the AGORA genus. -->


## Step 0: Checking taxonomic coverage

First, run the workflow with the `--dry-run` flag to prevent snakemake from submitting any jobs.  This will run the logic to match up the taxa you provided with those present in the pre-formated AGORA model.  Genera not able to matched will be written to `unmatchable_genera.csv` in your output directory. In cases of naming mismatches, create a csv file containing your and AGORA's genus names as columns `A` and `B`, and supply in the next step.

```
snakemake --config stage=tradeoff  abundances=$PWD/results/annotate/batch1_abundances_renamed.csv agora=$PWD/agora201_refseq216_species_1.qza medium=$PWD/western_diet_gut_agora.qza taxrank=species             cutoff=.0001 --directory $PWD/results/  --dry-run
```

## Step 1: Identifying an optimal tradeoff

We will utilize Snakemake config arguments to first submit jobs calculating the optimal tradeoff parameter for each sample.  The config can be modified to set the appropriate taxa abundance cutoff.

```
snakemake --config stage=tradeoff  abundances=$PWD/results/annotate/batch1_abundances_renamed.csv agora=$PWD/agora201_refseq216_species_1.qza medium=$PWD/western_diet_gut_agora.qza taxrank=species             cutoff=.0001 --directory $PWD/results/
```


## Step 2: Grow
Having decided on the tradeoff value based on the visualizations in `results/tradeoff.html`, you can now submit the growth jobs! This is done by changing the `stage` config argument, setting the `tradeoff` argument, and re-using the same output `--directory`.

```
snakemake --config stage=grow tradeoff=0.6  abundances=$PWD/results/annotate/batch1_abundances_renamed.csv agora=$PWD/agora201_refseq216_species_1.qza medium=$PWD/western_diet_gut_agora.qza taxrank=species             cutoff=.0001 --directory $PWD/results/
```

## Step 3: Analysis

From here, you can follow along the tutorial! Run the following in place of the `growth = grow(com, "models", medium, tradeoff=0.8, threads=2)` line:

```python
from micom.workflows import load_results
growth = load_results("results/growth.zip")
```
