ASSEMBLY="/work/rodrigues_2023/data/genome_assembly/MR_Cli"
ANNOTATION="02_rast_annotation_by_marinelle"

mkdir -p "$ASSEMBLY/$ANNOTATION/02_bowtie2_index"
cd "$ASSEMBLY/$ANNOTATION/02_bowtie2_index"
bowtie2-build "$ASSEMBLY/$ANNOTATION/all_ATEC_annotated_contigs.fa" MR_Cli_clc_rast