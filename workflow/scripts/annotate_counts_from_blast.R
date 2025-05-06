loge <- file(snakemake@log[["e"]], open = "wt")
logo <- file(snakemake@log[["o"]], open = "wt")
sink(logo, type = "output")
sink(loge, type = "message")
library(tidyverse)

print("Parsing annotation file")
counts <- read.csv(snakemake@input[["asv_counts"]])
for (col in c("sample_id", "asv", "count")){
    if (! col %in% colnames(counts)){
        stop("Counts file must contain columns named sample_id, asv, and count")
    }
}


blast_names <- c("ASVId", "taxid", "accession", "species", "query_length", "align_length", "nident", "pident", "bitscore", "evalue", "score")
results <- read.csv(snakemake@input[["blast"]], sep = "\t", header = FALSE, stringsAsFactors = FALSE, col.names = blast_names)

## following identifies the accessions with the best match, including tie scores
# find the tie score hits
# for the species, replace the first space with underscore, delete second space and beyond, then alphabetical order
blasted <- results %>%
  group_by(ASVId) %>%
  filter(bitscore == max(bitscore)) %>%
  arrange(accession) %>%
  # identify unique names, put them in alphabetical order, and merge them, and also add the % identity
  mutate(
    name = paste(paste(unique(accession), collapse = ";"), align_length[1], pident[1], sep = ";")
  ) %>%
  slice(1) %>%
  ungroup()



# this next section selects ASVs that should be added back if removed during chimera selection
# see how much of the ASV sequence BLAST was able to align
blasted$length_ratio <- blasted$query_length / blasted$align_length

# examine ASVs with a pident score 97% or better
blasted$pident_97 <- blasted$pident >= 97

# blasted_97<-blasted[blasted$pident_97==TRUE,]
# plot the results
# hist(blasted_97$length_ratio)
# it's a bimodal distribution, so let's use a cut-off of 1.1

print(str(blasted))

blasted$length_ratio_1.1 <- blasted$length_ratio <= 1.1
blasted$passed <- blasted$length_ratio_1.1 == TRUE & blasted$pident_97 == TRUE
print(table(blasted$passed))
# so these are the ASVs that blasted well enough to avoid chimera checking

# The output of blast is send to this file to be load in database.
write.table(blasted, snakemake@output[["detailed"]], sep = "\t", row.names = FALSE)



tmp <- counts %>%
    left_join(blasted %>% select(asv=ASVId, species=accession)) %>%
    group_by(sample_id) %>%
    mutate(abundance = count/sum(count)) %>%
    ungroup() %>%
    write.table(snakemake@output[["asv_abundances"]], row.names = FALSE, quote=FALSE, sep=",")
