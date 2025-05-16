stage=$1

case $stage in
    1|format16S)
	snakemake --config stage=format16S  agora_genomes=$PWD/.tests/genomes/ dbname="testdb"  --directory $PWD/mydb/
	;;
    2|annotate16S)
        snakemake --config stage=annotate16S  id=batch1 asv_counts=$PWD/.tests/asv_counts.csv asv_fasta=$PWD/.tests/asvs.fasta  model_manifest=$PWD/agora2_refseq_species.tsv db_path_nsq=$PWD/mydb/db/testdb.nsq  --directory $PWD/tmpresults/
	;;
    3|fixmodelid)
	echo "Dont be lazy, just make this a pandas script"
	singularity exec -B $PWD  docker://rocker/tidyverse:4.1 Rscript $PWD/workflow/scripts/01-genomeid_to_modelid_amplicon.R $PWD/tmpresults/annotate/batch1_abundances.csv $PWD/agora2_refseq_species.tsv $PWD/tmpresults/annotate/batch1_abundances_renamed.csv
	;;
    4|checkcoverage)
       python workflow/scripts/02-check_model_coverage.py --abundances $PWD/tmpresults/annotate/batch1_abundances_renamed.csv --agora ./agora201_refseq216_species_1.qza --outdir tmpresults/annotate/QC/
       ;;
    5|tradeoff)
	snakemake --config stage=tradeoff  abundances=$PWD/tmpresults/annotate/QC/postQC_abundances.csv  agora=$PWD/agora201_refseq216_species_1.qza medium=$PWD/western_diet_gut_agora.qza taxrank=species             cutoff=.0001 --directory $PWD/tmpresults/
	;;
    6|grow)
	snakemake --config stage=grow tradeoff=0.6  abundances=$PWD/tmpresults/annotate/batch1_abundances_renamed.csv agora=$PWD/agora201_refseq216_species_1.qza medium=$PWD/western_diet_gut_agora.qza taxrank=species             cutoff=.0001 --directory $PWD/tmpresults/
	 ;;
     * )
	 echo "not an option; select either format16S, annotate16S, fixmodelid, checkcoverage, tradeoff, or grow"
	 ;;
esac
