import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import pysam
import glob
import tqdm
import pybedtools
import argparse

'''
Author: Pratibha Jagannatha (pjaganna@ucsd.edu)
Input files: sorted and aligned long-reads (BAM) file, output directory, GTF file, directory containing NanoPlot output, QC cutoff
Output file: filtered reads, sorted and indexed

Example command:
python3 /path/to/LR-Ribo-STAMP/scripts/filter_bam_v2.py --out_dir /path/to/LR-Ribo-STAMP/test/input_dir/ --bam_path /path/to/LR-Ribo-STAMP/test/input_dir/Ribo-STAMP_HEK293T_rep1_chr1.sorted.bam --gtf /path/to/LR-Ribo-STAMP/test/ref/gencode.v19.annotation.chr1.gtf --qc_dir /path/to/LR-Ribo-STAMP/test/output_dir/qc/Ribo-STAMP_HEK293T_rep1_chr1/ --qc_cutoff 20

'''

def QC_read_filter(qc_directory, qc_score_cutoff):
    '''Get reads that need to be filtered out based on QC score from Nanoplot'''
    
    qc_df = pd.read_csv(os.path.join(qc_directory,'NanoPlot-data.tsv.gz'), sep='\t')
    qc_reads_to_remove = qc_df[qc_df['quals']<=qc_score_cutoff]['readIDs'].tolist()
    
    return list(set(qc_reads_to_remove))

def assign_strand(df):
    '''Function to define the mapping strand of a read'''
    
    if df["strand"] ==True:
        return "-"
    elif df["strand"] ==False:
        return "+"

def strand_read_filter(gtf_file, bam_file):
    '''Get reads that need to be filtered out based on mapping strand versus gene strand'''
    # get GTF BED information
    gtf_df = pd.read_csv(gtf_file, sep="\t", comment="#", header=None)
    gtf_df = gtf_df[gtf_df[2] == "gene"]  # Select rows with feature type "gene"
    gtf_df["name"] = gtf_df[8].str.extract(r'gene_name "([^"]+)"', expand=False)
    gtf_df = gtf_df.iloc[:, [0, 3, 4, 6,9]]
    gtf_df.columns = ["chrom", "start", "end", "strand", "name"]
    gtf_df["score"] = 0
    gtf_df = gtf_df[["chrom", "start", "end", "name", "score", "strand"]]
    # get BAM BED information 
    bam_df = pd.DataFrame()
    bam = pysam.AlignmentFile(bam_file)
    bam_df = pd.DataFrame([(read.reference_name, read.reference_start, read.reference_end, read.query_name, 0, read.is_reverse) for read in bam.fetch(until_eof=True)], columns=["chrom", "start", "end", "name", "score", "strand"])
    bam_df["strand"] = bam_df.apply(assign_strand, axis=1)
    
    # Convert the DataFrames to BedTool objects for overlap
    bam_bed = pybedtools.BedTool.from_dataframe(bam_df)
    gtf_bed = pybedtools.BedTool.from_dataframe(gtf_df)
    # Overlap coordinates
    intersected = bam_bed.intersect(gtf_bed, wo=True)
    intersected_df = intersected.to_dataframe(names=[]).reset_index()
    intersected_df.columns = ["read_chrom", "read_start", "read_end", "read_name","read_score", "read_strand", "gene_chrom", "gene_start", "gene_end", "gene_name",
           "gene_score", "gene_strand", "overlap"]
    # keeping the read name, gene combination with highest amount of overlap before removing strand mismatches
    intersected_df = intersected_df.sort_values(by=['read_name', 'overlap'], ascending=[True, False]).drop_duplicates(subset=['read_name'])
    # getting read names for ones with mismatched strand
    reads_strand_mismatch = intersected_df[intersected_df["gene_strand"] != intersected_df ['read_strand']]['read_name'].tolist()
    
    return list(set(reads_strand_mismatch))


def supplementary_read_filter(bam_file):
    '''Get names of supplementary alignment reads'''
    supplementary_readname = []
    samfile = pysam.AlignmentFile(bam_file, "rb")
    for read in samfile.fetch():
        if read.flag == 2048: 
            supplementary_readname.append(read.query_name)
        elif read.flag == 2064: 
            supplementary_readname.append(read.query_name)

    return list(set(supplementary_readname))


def filter_read(output_bam, bam_file, reads_to_filter):
    ''' '''
    key=bam_file.split('/')[-1].split('.sorted.bam')[0]
    bam = pysam.AlignmentFile(bam_file)
    obam = pysam.AlignmentFile(os.path.join(output_bam,key + '.readfiltered.bam'), "wb", template=bam)
    for q in bam.fetch(until_eof=True):
        if q.query_name not in reads_to_filter:
            if q.flag not in filter_out_flags:
                obam.write(q)
    obam.close()
    bam.close()
    pysam.sort("-o", os.path.join(output_bam,key + '.readfiltered.sorted.bam'), os.path.join(output_bam,key + '.readfiltered.bam'))
    pysam.index(os.path.join(output_bam,key + '.readfiltered.sorted.bam'))

parser = argparse.ArgumentParser(description='Identify read level edits. Make sure all split BAMs have index files prior to running script.')
parser.add_argument("--out_dir", type=str, default='./', help='directory where filtered BAM will be saved', required=True)
parser.add_argument("--bam_path", type=str, default='./', help='full path to unfiltered BAM', required=True)
parser.add_argument("--gtf", type=str, help='full path to GTF annotation (e.g. gencode.v19.annotation.gtf)', required=True)
parser.add_argument("--qc_dir", type=str, help='full path to nanoplot output.', required=True)
parser.add_argument("--qc_cutoff", type=int, default=20, help='cutoff threshold for Qscore', required=False)


args = parser.parse_args()
output_dir = args.out_dir
bam_file_path = args.bam_path
gtf_file_path = args.gtf
qc_dir = args.qc_dir
qc_threshold = args.qc_cutoff


filter_out_flags=[4,  # unmapped read
                  20, # unmapped read - reverse
                  256,# secondary alignment
                  272, # secondary alignment - reverse
                  2048, # supplementary alignment
                  2064 # supplementary alignment - reverse
                 ]

QC_reads = QC_read_filter(qc_dir, qc_threshold)
strand_reads = strand_read_filter(gtf_file_path, bam_file_path)
supp_reads = supplementary_read_filter(bam_file_path)

filter_out_reads = list(set(QC_reads)) + list(set(strand_reads)) + list(set(supp_reads))

print(f" Number of reads filtered based on Q score: {len(list(set(QC_reads)))}")
print(f" Number of reads filtered based on mismatched strand: {len(list(set(strand_reads)))}")
print(f" Number of reads filtered based on supplementary read flag:{len(list(set(supp_reads)))}")
print(f" Total reads filtered out: {len(list(set(filter_out_reads)))}")

filter_read(output_dir, bam_file_path, filter_out_reads)