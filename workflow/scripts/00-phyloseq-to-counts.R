library(dplyr)
args = commandArgs(trailingOnly=TRUE)
tmp <- readRDS(args[1])
tmp |> speedyseq::psmelt() |>
    select(asv=OTU, sample_id = Sample, count=Abundance) |>
    write.csv(args[2], row.names=FALSE, quote=FALSE)
