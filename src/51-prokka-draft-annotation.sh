ASSEMBLY="/work/rodrigues_2023/data/genome_assembly/MR_Cli/spades_run_03"
OUTDIR="$ASSEMBLY/02_prokka_annotation_usegenus"

prokka --compliant --centre UTSW --outdir $OUTDIR \
       --genus Escherichia --species coli --strain PEc --locustag PEc       \
       --evalue 0.001 --cpus 0 $ASSEMBLY/scaffolds.fasta