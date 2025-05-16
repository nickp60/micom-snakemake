rule all:
    input:
        f"annotate/{config['id']}_asv_blast_out.tsv",
        f"annotate/{config['id']}_abundances.csv",


rule legacy_blast_annotate:
    input:
        asv_fasta=config["asv_fasta"],
        db=config["db_path_nsq"],
    output:
        blast=f"annotate/{config['id']}_asv_blast_out.tsv",
    resources:
        mem_mb=32 * 1024,
    threads: 32
    params:
        dbdir = lambda wildcards, input: os.path.dirname(input.db),
        dbname = lambda wildcards, input: os.path.basename(os.path.splitext(input.db)[0])
    container:
        "docker://ghcr.io/vdblab/micro16s-blast:2.13.0a"
    shell:
        """
        export BLASTDB={params.dbdir}
        blastn \
            -num_threads {threads} \
            -db {params.dbname} \
            -max_target_seqs 500 \
            -query {input.asv_fasta} \
            -outfmt \"6 qseqid staxids saccver stitle qlen length nident pident bitscore evalue score\" \
            -out {output.blast}
        """




rule annotate_counts:
    input:
        asv_counts=config["asv_counts"],
        blast=f"annotate/{config['id']}_asv_blast_out.tsv",
    output:
        asv_abundances=f"annotate/{config['id']}_abundances.csv",
        detailed=f"annotate/{config['id']}_blast_passed_not_passed.txt",
    threads: 1
    resources:
        mem_mb=32 * 1024,
    container:
        "docker://rocker/tidyverse:4.1"
    log:
        e=f"logs/parse_blast_hits.e",
        o=f"logs/parse_blast_hits.o",
    script:
        "../scripts/annotate_counts_from_blast.R"
