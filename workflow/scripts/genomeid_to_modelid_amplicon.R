# This is an example script that assumes a lot!  YMMV
library(tidyverse)
args = commandArgs(trailing=TRUE)
if(length(args) != 3) stop("must have exactly 3 args : input.csv, agora2_refseq_species.tsv, and the path to the output.csv file name")

#if( !file.exists("agora2_refseq_species.tsv")) download.file("https://raw.githubusercontent.com/micom-dev/databases/refs/heads/main/recipes/agora2/manifests/agora2_refseq_species.tsv", "agora2_refseq_species.tsv")
key = read.csv(args[2], sep="\t")

# here we fix the one that didn't match on species >:(
key$species <- ifelse(grepl("lavalensis", key$strain), key$strain, key$species)

tmp <- read.csv(args[1])

# merge with manifest
tmp3 <- left_join(tmp %>% rename(id = species),
                  key %>% select(id, species) %>% distinct())
print("tmp3")
print(head(tmp3))
# In summary,
# - we used sylph with the AGORA genomes to get a taxonomic profile
# - we merged that with the agora2_refseq_species.tsv to get the full taxonomy; we merge on the species after stripping out brackets of the key
# - (ignore!) we then replace special char with underscores to match the file names present in agora201_refseq216_species_1.qza
 write.table(tmp3 %>%
             select(sample_id, abundance, species)%>%
             filter(abundance > 0) %>%
             group_by(sample_id) %>%
             arrange(sample_id, desc(abundance)) %>%
             ungroup(),
             args[3], row.names=FALSE, sep=",", quote=FALSE)
