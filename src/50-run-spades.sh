#!/bin/bash
# This script runs SPAdes using a YAML configuration file as input

# Set the input and output directory paths
#OUTPUT_DIR="/work/rodrigues_2023/data/genome_assembly/20230503/01_spades_output"
OUTPUT_DIR="/work/rodrigues_2023/data/genome_assembly/MR_Cli/spades_run_03"

# Set the path to the SPAdes executable
SPADES_PATH="/work/rodrigues_2023/tools/SPAdes-3.15.5-Linux/bin/spades.py"

# Set the path to the YAML configuration file
CONFIG_FILE="/work/rodrigues_2023/data/genome_assembly/20230503/spades_input.yml"

# Run SPAdes with the specified parameters using the YAML configuration file as input
#$SPADES_PATH -t 16 -m 370 -o $OUTPUT_DIR -k 21,33,55,77,99,127 --dataset $CONFIG_FILE --isolate
$SPADES_PATH -t 8 -m 250 -o $OUTPUT_DIR -k 21,33,55,77,99,127 --isolate     \
--pe1-1 "/work/rodrigues_2023/data/raw/isolates/MR_Cli1_1.fq.gz" \
--pe1-2 "/work/rodrigues_2023/data/raw/isolates/MR_Cli1_2.fq.gz"
# 