#!/usr/bin/env python

import os
import pandas as pd
import pysam
import numpy as np
import pickle
from collections import defaultdict
import seaborn as sns
from tqdm import tqdm, trange
from math import ceil
import argparse

"""
Author: Pratibha Jagannatha (pjaganna@ucsd.edu)
Input files: filtered sorted and aligned long-reads (BAM) file, output directory, directory with isoquant output, sample name
Output file: split BAM files, dictionary in pickle format to map genes and isoforms back to splits

Script to split the filtered BAM file into smaller chunks for faster processing.
"""

def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Split BAM file based on isoform read assignments.")
    parser.add_argument("--isoquant_dir", type=str, help="Directory containing isoform quantification (Isoquant) results")
    parser.add_argument("--genome_bam", type=str, help="Path to the genome BAM file")
    parser.add_argument("--output_dir", type=str, help="Output directory for split BAM files")
    parser.add_argument("--sample_name", type=str, help="Name of the sample")
    return parser.parse_args()

def load_isoquant_table(isoquant_dir, sample_name):
    """Load the isoquant table and filter unique isoforms."""
    isoquant_table = pd.read_csv(os.path.join(isoquant_dir, f'{sample_name}.read_assignments.tsv'), sep='\t', comment='#', header=None)
    isoquant_table_unique = isoquant_table[(isoquant_table[5]=='unique') | (isoquant_table[5]=='unique_minor_difference')] # only use uniquely assigned reads
    transcript_to_gene = dict(zip(isoquant_table_unique.drop_duplicates(subset=3)[3], isoquant_table_unique.drop_duplicates(subset=3)[4]))
    iso_read_assignment = isoquant_table_unique.groupby(3)[0].apply(list).reset_index()
    iso_read_assignment = iso_read_assignment.rename(columns={3 :'isoform', 0:'read_id'})
    iso_read_assignment = iso_read_assignment[iso_read_assignment['isoform'].str.contains('ENST')]
    iso_read_assignment['gene'] = iso_read_assignment['isoform'].map(transcript_to_gene)
    iso_read_assignment['num_reads'] = iso_read_assignment.apply(lambda x: len(x['read_id']), axis=1)
    print(f'Reads assigned to {len(set(iso_read_assignment["isoform"].tolist()))} transcripts')
    iso_read_assignment['iso_group'] = iso_read_assignment.sort_values(['gene']).groupby('gene').cumcount() + 1
    return iso_read_assignment

def define_splits(isoform_read_info):
    """Define the split groups for isoforms."""
    S = max(isoform_read_info['iso_group'])
    isoform_read_info['iso_group'] = isoform_read_info['iso_group'].apply(str)

    split_s_dict = {}
    for i in range(1, int(S) + 1):
        s = isoform_read_info[isoform_read_info['iso_group'] == str(i)]
        split_s = np.array_split(s, ceil(s.shape[0] / 100))
        count = 0
        for l in range(0, len(split_s)):
            key = 'split_bam_{}_{}'.format(str(i), str(l))
            split_s_dict[key] = split_s[l]
            count += len(split_s[l])

        if count != int(isoform_read_info['iso_group'].value_counts().loc[str(i)]):
            print('Some iso info is missing')
            break
    return split_s_dict

def build_split_dictionary(split_s_dict):
    """Build a dictionary of split_bam keys : set of read names as values."""
    splits = defaultdict(set)
    for key in tqdm(split_s_dict.keys()):
        reads_df = pd.DataFrame([read for sublist in split_s_dict[key].loc[:, 'read_id'].tolist() for read in sublist])
        if reads_df.shape[0] >= 1:
            splits[key] = set(reads_df[0])
        else:
            print('Split {} has no reads.'.format(str(b)))
    return splits

def build_split_read_object_dictionary(split_readname_dict, possorted_bam):
    """Build a dictionary of split_bam keys : list of read objects as values."""
    mapped_reads = int(pysam.view("-c", possorted_bam).strip('\n'))
    progress = trange(mapped_reads, position=0, leave=True)
    records_dict = defaultdict(list)
    bam = pysam.AlignmentFile(possorted_bam, "rb")
    for q in bam.fetch(until_eof=True):
        for split_key, qnames in split_readname_dict.items():
            if q.query_name in qnames:
                records_dict[split_key].append(q)
                qnames.remove(q.query_name)
        progress.update(1)
    return records_dict

def create_split_bam_files(output_bam_split, records_dict, possorted_bam):
    """Create split BAM files."""
    bam = pysam.AlignmentFile(possorted_bam, "rb")
    progress = trange(len(records_dict.keys()))
    for key, reads in records_dict.items():
        if len(reads) >= 1:
            obam = pysam.AlignmentFile(os.path.join(output_bam_split, key + '.bam'), "wb", template=bam)
            for read in reads:
                obam.write(read)
            obam.close()
            if int(pysam.view("-c", os.path.join(output_bam_split, key + '.bam')).strip('/n')) != len(reads):
                print(f"{key} | reads in BAM do not match reads assigned to split for this sample")
            pysam.sort("-o", os.path.join(output_bam_split, key + '.sorted.bam'), os.path.join(output_bam_split, key + '.bam'))
            pysam.index(os.path.join(output_bam_split, key + '.sorted.bam'))
        else:
            print(f"{key} | has no reads.")
        progress.update(1)

if __name__ == "__main__":
    args = parse_arguments()

    iso_read_assignment = load_isoquant_table(args.isoquant_dir, args.sample_name)
    iso_read_assignment.to_csv(os.path.join(args.output_dir, f'{args.sample_name}_transcript_counts.txt'), sep='\t', index=False)
    split_dictionary = define_splits(iso_read_assignment)

    # Output pickle file
    filename = os.path.join(args.output_dir, f"split_bam_{args.sample_name}")
    with open(filename, 'wb') as outfile:
        pickle.dump(split_dictionary, outfile)

    split_to_readname_dict = build_split_dictionary(split_dictionary)
    records_to_create_splits = build_split_read_object_dictionary(split_to_readname_dict, args.genome_bam)
    create_split_bam_files(args.output_dir, records_to_create_splits, args.genome_bam)
    
    with open(os.path.join(args.output_dir, 'completed.log'), 'w') as fp:
        pass
