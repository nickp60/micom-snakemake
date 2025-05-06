# micom-snakemake
A snakemake workflow for running MICOM  https://micom-dev.github.io/micom/

## 16S amplicon data

After your sequences-by-counts table has been generated (DADA2, DEBLUR), extract the ASV/OTU sequences as a fasta file.  These can be annotated by blasting against a database created from the AGORA genomes.

### Download AGORA2 genomes
Download AGORA2 genomes as described [here](https://github.com/vdblab/resources/tree/main/workflow/sylphflux)

### extract 16S and build BLAST database




## Shotgun metagenomic data
### Get a Shotgun Taxonomic Profile
We use sylph to get taxonomic profiles against a database constructed from the AGORA2 genomes. See the workflow [here](https://github.com/vdblab/resources/tree/main/workflow/sylphflux) for downloading the genomes and formatting into a sylph sketch.  Extract a table from the results that is in "tall" format, with columns for `sample_id`, `species/genus`, and `abundance`.

###
