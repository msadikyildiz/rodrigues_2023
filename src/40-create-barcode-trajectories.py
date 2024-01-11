import pysam
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse

parser = argparse.ArgumentParser(description='Finds barcoded reads from input bam file.')
parser.add_argument('--input-dir', help='Input directory path where barcode counts live.', required=True)
parser.add_argument('--metadata', help='Tab seperated metadata file where sample names and timeline provided.', required=True)
parser.add_argument('--output-dir', help='Output folder where you want trajectories to be printed.', required=True)

args = vars(parser.parse_args())

meta = pd.read_csv(args['metadata'], sep='\t', index_col=0)

for sample in meta.columns:
    for day in meta.index:
