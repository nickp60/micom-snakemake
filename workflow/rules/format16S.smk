import yaml
import sys
import os
import glob


envvars:
    "TMPDIR",




onstart:
    with open("annotate_config_used.yaml", "w") as f:
        yaml.dump(config, f)


localrules:
    all,


# from https://stackoverflow.com/questions/3703276
def is_gz_file(filepath):
    with open(filepath, 'rb') as test_f:
        return test_f.read(2) == b'\x1f\x8b'


# see note in resources repo and this repo's REDADME; empty files exist due to ncbi's 'suppressed' status of certain accessions
# the workflow sometimes donwloads html for malformed links; this checks for gzipped files
files = glob.glob(config["agora_genomes"] + "*.fna.gz")
print(len(files))
files = [x for x in files if os.path.getsize(x) > 0]
print(len(files))
files = {os.path.basename(x).replace(".fna.gz", ""):x for x in files if is_gz_file(x)}
print(len(files.keys()))

rule all:
    input:
        multiext(f"db/{config["dbname"]}", ".nsq", ".nin",  ".nhr")


# rule barnapp_agora:
#     params:
#         dbdir = config["agora_genomes"],

#     output:
#         fasta = "db/agora16S.fa",
#       #  indfnas = temp(directory("barrnaptmp")),
#     container: "docker://chrishah/barrnap:0.9",
#     threads: 8
#     shell:"""
#     mkdir -p barrnaptmp/
#     # see note in resources repo; empty files exist due to ncbi's 'suppressed' status of certain accessions
#     find {params.dbdir} -name "*.fna.gz" -size +0 -print  | while read fnagz;
#     do
#     # the workflow sometimes donwloads html for malformed links; this checks for gzipped files
#     if gunzip --test "${{fnagz}}" 2>/dev/null 1>/dev/null; then
#     basen=$(basename $fnagz | sed "s|.fna.gz||g")
#     if [ ! -f "barrnaptmp/${{basen}}.fa" ]; then
#     zcat $fnagz | barrnap -o barrnaptmp/${{basen}}.fa --threads {threads} -  > /dev/null
#     fi
#     fi
#     done
#     cat barrnaptmp/*fa > {output.fasta}
#     """
def get_genome_fna(wc):
    return files[wc.name]

rule barnapp_agora:
    input:
        fasta = get_genome_fna,
    output:
        fasta = "barrnaptmp/{name}.fna"
    container: "docker://chrishah/barrnap:0.9",
    threads: 1
    shell:"""
    zcat {input.fasta} | barrnap -o {output.fasta} --threads {threads} -  > /dev/null
    """

rule agg_and_reheader:
    """ we rename the headers so we get the full species name in the blast output
    """
    input:
        expand("barrnaptmp/{name}.fna", name=files.keys())
    output:
        fasta = f"db/{config['dbname']}.fa"
    container: "docker://ghcr.io/vdblab/seqkit:2.10.0a"
    # fbne and nr are special seqkit syntax to represent the extensionless file name and line number, respectively
    shell: """
    for x in {input};
    do
    seqkit replace $x -p '^' -r '{{fbne}} {{nr}} ' >> {output.fasta}
    done
    """

rule agora_blast_db:
    container:
        "docker://ncbi/blast:2.7.1"
    input:
        fasta = f"db/{config['dbname']}.fa"
    output:
        multiext(f"db/{config['dbname']}", ".nsq", ".nin",  ".nhr")
    resources:
        mem_mb=8 * 1024,
        runtime=4*60,
    params:
        dbpath = lambda wildcards, input: input.fasta.replace(".fa", ""),
        dbname = lambda wildcards, input: os.path.basename(input.fasta.replace(".fa", "")),
    shell:"""
    makeblastdb  -in {input.fasta} -dbtype nucl -out {params.dbpath} -title {params.dbname}
    find db
    """
