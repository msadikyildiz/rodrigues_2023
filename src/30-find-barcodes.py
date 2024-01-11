import pysam
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse

def find_base_in_alignment(alignment: pysam.AlignedSegment, pos: int):
    idx_q = 0
    idx_r = pos - alignment.reference_start
    
    seq = alignment.query_sequence

    # if alignment.is_reverse:
    #     seq = alignment.get_forward_sequence()
    # else:
    #     seq = alignment.query_sequence
    
    if seq is None:
        return None
    
    for op, l in alignment.cigartuples:
        ref_consumed = op in {0, 2, 3, 7, 8}
        query_consumed = op in {0, 1, 4, 7, 8}
        
        if ref_consumed:
            idx_r -= l
        if query_consumed:
            idx_q += l
        
        if idx_r < 0:
            if query_consumed:
                # base is in query between idx_q-l , idx_q
                base = seq[idx_q + idx_r - 1]
                return base
            else:
                # position has been deleted
                return None


parser = argparse.ArgumentParser(description='Finds barcoded reads from input bam file.')
parser.add_argument('--bam', help='Input bam file path', required=True)
parser.add_argument('--fasta', help='Input reference fasta file path', required=True)
parser.add_argument('--contig', help='Input reference contig', required=True)
parser.add_argument('--barcode-start', help='Barcode start position', required=True)
parser.add_argument('--barcode-stop', help='Barcode stop position', required=True)
parser.add_argument('--maximum-missing-coverage', help='How many N positions is allowed to be missing?', required=True)
parser.add_argument('--output', help='Specify the output folder where counts will be printed.', required=True)

args = vars(parser.parse_args())

bam = pysam.AlignmentFile(args["bam"], "rb")
ref = pysam.FastaFile(args["fasta"])

#barcode_positions = [134,135,136,139,140,141,142,145,146,147,148]
#barcode1 = [133,147]
#barcode2 = [152,166]
barcode_start = int(args["barcode_start"])
barcode_stop = int(args["barcode_stop"])
barcode_print_stop = 171
# Find N positions within barcode region
refseq = ref.fetch(args["contig"], start=barcode_start, end=barcode_stop)
barcode_var_positions = barcode_start + np.where(pd.Series(list(refseq))=='N')[0]
print("Variable positions on the specified barcode region:", barcode_var_positions)

barcoded_reads=[]
for read in bam.fetch('amplicon_bgiseq', start=barcode_start, stop=barcode_stop):
    if ((read.reference_start + read.query_alignment_start) < barcode_start) and \
        ((read.reference_start + read.query_alignment_end) >= barcode_stop):

        barcoded_read={'reference_start': read.reference_start,
                       'is_read1': 1 if read.is_read1 else 2,
                       'strand': '+' if read.is_reverse else '-',
                       'query_alignment_start': read.query_alignment_start,
                       'query_alignment_end': read.query_alignment_end}
        full_barcode=[]
        var_barcode=[]
        for i in range(barcode_start, barcode_print_stop):
            base = find_base_in_alignment(read,i)
            if base is None:
                base = '_'   
            full_barcode.append(base)
            if i in barcode_var_positions:
                var_barcode.append(base)

        barcoded_read['full_barcode'] = "".join(full_barcode)
        barcoded_read['variable_barcode'] = "".join(var_barcode)
        barcoded_reads.append(barcoded_read)

# print out barcoded reads
df = pd.DataFrame(barcoded_reads)
df.to_csv(args["output"]+'.tsv', sep='\t', index=False)

# count barcodes
counts={}
for bar in df['variable_barcode']:
    if bar.count('_') <= int(args['maximum_missing_coverage']):
        if bar in counts:
            counts[bar] += 1
        else:
            counts[bar] = 1

# print barcode counts
with open(args["output"]+".counts", "w") as f:
    for bar in sorted(list(counts.items()), key=lambda x: x[1], reverse=True):
        f.write(f"{bar[0]}\t{bar[1]}"+"\n")