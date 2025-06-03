# This is an example script that assumes a lot!  YMMV
library(tidyverse)
args = commandArgs(trailing=TRUE)
if(length(args) != 3) stop("must have exactly 3 args : wide_sylph_output.csv, agora_manifest.tsv, and output.csv file name")

key = read.csv(args[2], sep="\t")
# here we fix the one that didn't match on species >:(
key$species <- ifelse(grepl("lavalensis", key$strain), key$strain, key$species)

tmp <- read.csv(args[1]) %>% filter(startsWith(clade_name, "d__Bacteria")) %>%
select(-taxids) %>% pivot_longer(cols=-clade_name, names_to="sample_id", values_to="abundance") %>%
mutate(species = gsub(".*s__(.*)\\|t__.*", "\\1", clade_name)) %>%
    group_by(sample_id, species) %>%
    summarize(abundance=sum(abundance)) %>%
group_by(sample_id) %>%
mutate(relabundance=abundance/sum(abundance)) %>%
ungroup() %>%
mutate(sample_id = gsub("^X", "", sample_id)) %>%
       select(sample_id, species, abundance=relabundance)
print(str(tmp))
# merge with manifest
tmp3 <- left_join(tmp %>% rename(species_clean = species),
                  key %>% mutate(species_clean = gsub("\\[|\\]", "", species)) %>% select(species_clean, species) %>% distinct())
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
