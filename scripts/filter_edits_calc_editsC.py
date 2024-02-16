import os
import pybedtools
import pandas as pd
import numpy as np
import glob
import tqdm
import seaborn as sns
import matplotlib.pyplot as plt
import random
from Bio import SeqIO
from Bio.Seq import Seq
import pysam
import argparse


"""
Author: Pratibha Jagannatha (pjaganna@ucsd.edu)
Input files: tab-delimited sample file, SNP file, gene fasta, isoform fasta, output directory
Output file: filtered edits, EditsC tables

Script to filter edits and calculate EditsC.

Example: python3 /path/to/LR-Ribo-STAMP/scripts/filter_edits_calc_editsC.py --sample_file /path/to/sample_file.txt --snp_file /path/to/hg19.commonSNPs147.bed3 --gtf_file/path/to/gencode.v19.annotation.gtf --gene_coordinates_fasta /path/to/merged_exons_utrs.genes.hg19.fasta --transcript_coordinates_fasta /path/to/merged_exons_utrs.transcript.hg19.fasta --output_dir /path/to/output

"""

def get_sample_information(sample_info_path):
    '''parse sample file to get data information'''
    sample_info = pd.read_csv(sample_info_path, sep='\t')
    samples_dict = {}
    cs = list(set(sample_info['condition'].tolist()))
    for i in range(len(cs)):
        samples_dict[cs[i]] = dict(zip(sample_info[sample_info['condition']==cs[i]].rep, sample_info[sample_info['condition']==cs[i]].sample_name))
    sample_paths = dict(zip(sample_info.sample_name, sample_info.path))

    return samples_dict, sample_paths

def concat_splits_all_pos_edits(sample_dict):
    '''concatenate CDS edits called from each split BAM'''
    all_sample_pos_edits = {}
    for key in sample_dict:
        print(key)
        all_pos_edits_df = pd.DataFrame()
        empty_splits = []
        bed_files = [i for i in glob.glob(f'{map_sample_path[sample_dict[key]]}/output_dir/*_ct_*.bed') if '_transcript' not in i]
        for p in tqdm.tqdm(bed_files):
            if '_transcript' not in p:
                if os.path.getsize(p) != 0:
                    split_pos_edit_df = pd.read_csv(p, sep='\t', header=None)
                    split_pos_edit_df['split_label'] =  p.split('/')[-1].split('_ct_')[0]
                    split_pos_edit_df['key'] = split_pos_edit_df.apply(lambda x: f'{x[0]}:{str(x[1])}:{str(x[2])}', axis=1)
                    split_pos_edit_df['location'] = p.split('/')[-1].split('.bed')[0].split('_ct_')[-1]
                    split_pos_edit_df['conversion'] = split_pos_edit_df.apply(lambda x: int(x[4].split(',')[0]), axis=1)
                    split_pos_edit_df['coverage'] = split_pos_edit_df.apply(lambda x: int(x[4].split(',')[1]), axis=1)
                    split_pos_edit_df['edit_fraction'] = split_pos_edit_df.apply(lambda x: x['conversion']/x['coverage'], axis=1)
                    split_pos_edit_df = split_pos_edit_df[split_pos_edit_df['location']=='CDS']#### ONLY KEEPING CDS EDITS
                    all_pos_edits_df = pd.concat([all_pos_edits_df, split_pos_edit_df])
                else:
                    empty_splits.append(p.split('/')[-1])

        print(' '.join(empty_splits))
        all_sample_pos_edits[key] = all_pos_edits_df
        
    return all_sample_pos_edits
    

def get_common_sites(all_edits):
    '''get common sites across samples and replicates'''
    intersection_set = set()
    for sample in all_edits:
        for rep in all_edits[sample]:
            if intersection_set == set():
                intersection_set = set(all_edits[sample][rep]['key'].tolist())
            else:
                intersection_set = intersection_set.intersection(set(all_edits[sample][rep]['key'].tolist()))
    return intersection_set

def count_cytosines(fasta_file):
    '''Count cytosines at gene level'''
    cytosine_counts = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        gene_ids = record.id.split(":")[0].split(',')  # Adjust based on your FASTA header format
        sequence = str(record.seq).upper()
        for gene_id in gene_ids:
            if gene_id in cytosine_counts:
                cytosine_counts[gene_id] += sequence.count('C')
            else:
                cytosine_counts[gene_id] = sequence.count('C')
    return cytosine_counts

def count_cytosines_transcript(fasta_file):
    '''Count cytosines at isoform level'''
    cytosine_counts = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        ids = record.id.split(":")[0].split(',')  # Adjust based on your FASTA header format
        sequence = str(record.seq).upper()
        for transcript in ids:
            if transcript in cytosine_counts:
                cytosine_counts[transcript] += sequence.count('C')
            else:
                cytosine_counts[transcript] = sequence.count('C')
    return cytosine_counts


def calculate_editsC(keep_edits, c_counts):
    '''Calculate EditsC at gene level'''
    gene_editsC = {}
    for rep in ['rep1', 'rep2', 'rep3']:
        edit_count_gene = keep_edits[rep].drop_duplicates(subset=[0,1,2,5, 'gene_id']).groupby('gene_id').count()[[0]]
        edit_count_gene['num_c'] = edit_count_gene.index.map(c_counts)
        edit_count_gene = edit_count_gene.rename(columns={0:'edited_c'})
        edit_count_gene['editsC'] = edit_count_gene['edited_c']/edit_count_gene['num_c']
        gene_editsC[rep] = edit_count_gene
    return gene_editsC

def calculate_editsC_transcript(keep_edits, c_counts):
    '''Calculate EditsC at isoform level'''
    transcript_editsC = {}
    for rep in ['rep1', 'rep2', 'rep3']:
        edit_count_transcript = keep_edits[rep].groupby('transcript_id').count()[[0]]
        edit_count_transcript['num_c'] = edit_count_transcript.index.map(c_counts)
        edit_count_transcript = edit_count_transcript.rename(columns={0:'edited_c'})
        edit_count_transcript['editsC'] = edit_count_transcript['edited_c']/edit_count_transcript['num_c']
        transcript_editsC[rep] = edit_count_transcript
    return transcript_editsC




def parse_args():
    parser = argparse.ArgumentParser(description='Process input files and output directory.')
    parser.add_argument('--sample_file', type=str, required=True, help='Tab-delimited file containing sample name, replicate, conditions, and path to working directory')
    parser.add_argument('--snp_file', type=str, required=True,
                        help='Path to the SNP file. Example: /path/to/hg19.commonSNPs147.bed3')
    parser.add_argument('--gtf_file', type=str, required=True,
                        help='Path to the GTF file. Example: /path/to/gencode.v19.annotation.gtf')
    parser.add_argument('--gene_coordinates_fasta', type=str, required=True,
                        help='Path to the gene coordinates fasta file. Example: /path/to/merged_exons_utrs.genes.hg19.fasta')
    parser.add_argument('--transcript_coordinates_fasta', type=str, required=True,
                        help='Path to the transcript coordinates fasta file. Example: /path/to/merged_exons_utrs.transcript.hg19.fasta')
    parser.add_argument('--output_dir', type=str, required=True,
                        help='Path to the output directory. Example: /path/to/output/')

    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = parse_args()
    sample_info = args.sample_file
    snp_file = args.snp_file
    gtf_file =  args.gtf_file
    gene_coordinates_fasta = args.gene_coordinates_fasta
    transcript_coordinates_fasta = args.transcript_coordinates_fasta
    output_dir =  args.output_dir

    map_sample_condition, map_sample_path = get_sample_information(sample_info)
    all_pos_edits = {}
    for s in map_sample_condition:
        all_pos_edits[s] = concat_splits_all_pos_edits(map_sample_condition[s])
    common_sites = get_common_sites(all_pos_edits)

    # process SNP data
    snps = pd.read_csv(snp_file, sep='\t', header=None)
    snps['key'] = snps.apply(lambda x: f'{x[0]}:{str(x[1])}:{str(x[2])}', axis=1)
    snps['length'] = snps[2] - snps[1]
    snps_list = snps['key'].tolist()

    # remove common sites, SNPs
    all_pos_edits_common_removed = {key: {} for key in all_pos_edits.keys()}
    all_pos_edits_common_all_rps2_snp_removed = {key: {} for key in all_pos_edits.keys()}
    keep_edits = {key: {} for key in all_pos_edits.keys()}
    for s in all_pos_edits:
        for r in all_pos_edits[s]:
            all_pos_edits_common_removed[s][r] = all_pos_edits[s][r][~all_pos_edits[s][r]['key'].isin(common_sites)]
            all_pos_edits_common_all_rps2_snp_removed[s][r] = all_pos_edits_common_removed[s][r][~all_pos_edits_common_removed[s][r]['key'].isin(snps_list)]
            edits_df = all_pos_edits_common_all_rps2_snp_removed[s][r].copy()
            edits_df['gene_id'] = edits_df.apply(lambda x: x[3].split(':')[0], axis=1)
            edits_df['transcript_id'] = edits_df.apply(lambda x: x[3].split(':')[1], axis=1)
            keep_edits[s][r] = edits_df

    print('FILTERED EDITS:')
    for s in keep_edits:
        for r in keep_edits[s]:
            print(f'{s} {r} : {keep_edits[s][r].shape[0]}')
            keep_edits[s][r].to_csv(os.path.join(output_dir, f'{s}_{r}_filtered_edits.txt'), sep='\t', index=False)

    gene_cytosine_counts = count_cytosines(gene_coordinates_fasta)
    transcript_cytosine_counts = count_cytosines_transcript(transcript_coordinates_fasta)

    # Calculate EditsC
    gene_editsC = {}
    transcript_editsC = {}
    for s in keep_edits:
        gene_editsC[s] = calculate_editsC(keep_edits[s], gene_cytosine_counts)
        transcript_editsC[s] = calculate_editsC_transcript(keep_edits[s], transcript_cytosine_counts)

    # Output results
    for s in gene_editsC:
        for r in gene_editsC[s]:
            gene_editsC[s][r][['editsC']].dropna().to_csv(os.path.join(output_dir,f'{s}_{r}_gene_editsC.txt'), sep='\t')
            transcript_editsC[s][r][['editsC']].dropna().to_csv(os.path.join(output_dir, f'{s}_{r}_isoform_editsC.txt'), sep='\t')




























