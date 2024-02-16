#!/usr/bin/env python

import pandas as pd
import math
import pysam
import pybedtools
import os
from collections import defaultdict
import numpy as np
import pysam
from multiprocessing import  Pool
from multiprocessing import  cpu_count
import datetime
import pickle
from scipy.special import betainc
import argparse

'''
Author: Pratibha Jagannatha (pjaganna@ucsd.edu)
Input files: directory containing split BAM files (in_dir), output directory (generally same as in_dir), GTF, Fasta, name for split dictionary
Output file: edit files

This script is used to identify read-level edits when using annotated transcripts. It should be run following split_bam_isoquant.py as some of the output files from that script are used in this one.

Example of how to run this script:

python3 /path/to/read_level_quant_se_ct_annotated.py --in_dir /path/to/output_dir/ --out_dir /path/to/output_dir/ --gtf_path /path/to/refs/gencode.v19.annotation.gtf --ref /path/to/refs/hg19.fa --split_pickle /path/to/split_bam_dictionary --split split_bam_1_0

'''

def get_position_matrix(isoform_group_gff, stepper='nofilter', debug_read="", min_base_quality=0):
    '''This function was primarily authored by Brian Yee. It returns a dictionary of presence or absense of edits at a give position within a read.'''
    bam=sorted_bam
    infile = pysam.AlignmentFile(bam, "rb", reference_filename=reffile)
    edits_df = pd.DataFrame()
    pos_edits = pd.DataFrame()
    for row, index in isoform_group_gff.iterrows():
        chrom, start, end, strand_assign = isoform_group_gff.loc[row,0],isoform_group_gff.loc[row,3],isoform_group_gff.loc[row,4], isoform_group_gff.loc[row,6]
        print(chrom, start, end, strand_assign) ##########
        reference = pybedtools.BedTool.seq([chrom,0,end],reffile)
        count = start # running counter for each position added
        alphabet = {}
        edits = defaultdict(dict)
        positions = []
        offset = 0
        max_offset = 0
        MAX_DEPTH = 10000000
        for pileupcolumn in infile.pileup(chrom, start, end, stepper=stepper, max_depth=MAX_DEPTH, min_base_quality=min_base_quality):
            if pileupcolumn.pos >= start:  # pileupcolumn.pos is the lowest genomic coordinate (0-based) across all reads at <start>
                st = ""  # Store all bases foreach read at ONE position (count)
                if count >= (end) or pileupcolumn.pos >= end:  # I think this works because of sorted reads?
                    break
                """
                if there is no read coverage at the beginning positions
                """
                while(count < pileupcolumn.pos):
                    alphabet['A'] = 0
                    alphabet['T'] = 0
                    alphabet['C'] = 0
                    alphabet['G'] = 0
                    alphabet['del'] = 0
                    alphabet['ins'] = 0
                    alphabet['ref'] = reference[pileupcolumn.reference_pos].upper()
                    positions.append(alphabet)
                    alphabet = {}
                    count = count + 1
                
                """
                count all the bases at each position, store in string (st)
                """
                for pileupread in pileupcolumn.pileups: # for each pileup read
                    edits[pileupread.alignment.query_name][count] = 0
                    if (not pileupread.is_del) and (not pileupread.is_refskip):
                        """
                        THIS BLOCK DETERMINES IF C>T or G>A (or A>G or T>C)EDIT EXISTS
                        """
                        #print(pileupread.query_position) #####
                        if edit_type == 'ct':
                            if strand_assign == '+' and reference[count].upper() == 'C' and pileupread.alignment.query_sequence[pileupread.query_position].upper() == 'T':
                                edits[pileupread.alignment.query_name][count] = 1
                            elif strand_assign == '-' and reference[count].upper() == 'G' and pileupread.alignment.query_sequence[pileupread.query_position].upper() == 'A':
                                edits[pileupread.alignment.query_name][count] = 1

                        elif edit_type == 'ag':
                            if strand_assign == '+' and reference[count].upper() == 'A' and pileupread.alignment.query_sequence[pileupread.query_position].upper() == 'G':
                                edits[pileupread.alignment.query_name][count] = 1
                            elif strand_assign == '-' and reference[count].upper() == 'T' and pileupread.alignment.query_sequence[pileupread.query_position].upper() == 'C':
                                edits[pileupread.alignment.query_name][count] = 1

                        st = st + pileupread.alignment.query_sequence[pileupread.query_position].upper()

                        if(pileupread.alignment.query_name==debug_read):
                            print("st at refpos: {} and querypos: {}: {}".format(count, pileupread.query_position, st))
                    elif pileupread.is_del:
                        st = st + 'd'
                        if(pileupread.alignment.query_name==debug_read):
                            print("st at refpos: {} and querypos: {}: {}".format(count, pileupread.query_position, st))
                    elif pileupread.is_refskip:
                        st = st + 'i'
                        if(pileupread.alignment.query_name==debug_read):
                            print("st at refpos: {} and querypos: {}: {}".format(count, pileupread.query_position, st))
                    else:
                        print("THIS READ", pileupread.alignment.query_name, count, pileupcolumn.reference_pos, pileupread.query_position)
                        print( pileupread.is_del, pileupread.is_refskip, pileupread.indel)
                        st = st + '-'
                        if(pileupread.alignment.query_name==debug_read):
                            print("st at refpos: {} and querypos: {}: {}, [{}]".format(count, pileupread.query_position, st, pileupread.alignment.query_sequence[76]))
                            
                """
                count number of occurrences of each base
                """
                alphabet['A'] = st.count('A')
                alphabet['T'] = st.count('T')
                alphabet['C'] = st.count('C')
                alphabet['G'] = st.count('G')
                alphabet['del'] = st.count('d')
                alphabet['ins'] = st.count('i')
                alphabet['ref'] = reference[pileupcolumn.reference_pos].upper()
                count = count + 1
                # print(alphabet)

                positions.append(alphabet)
                alphabet = {}
        """
        If there are positions in the end without read coverage
        """
        while(count < end):
            # count = count + 1
            alphabet['A'] = 0
            alphabet['T'] = 0
            alphabet['C'] = 0
            alphabet['G'] = 0
            alphabet['del'] = 0
            alphabet['ins'] = 0
            alphabet['ref'] = reference[count].upper()

            for read_name in edits.keys():
                edits[read_name][count] = np.nan

            count = count + 1
            positions.append(alphabet)
            #pybedtools.helpers.cleanup()
        
        df = pd.DataFrame(positions)
        df.index = df.index + start  # re-align index to match genomic coordinates.
        df['chrom'] = chrom
        df['isoform'] = isoform_group_gff.loc[row,'isoform']
        df['gene'] = isoform_group_gff.loc[row,'gene']
        df['strand'] = strand_assign
        if edit_type == 'ct':
            if strand_assign == '+':
                df = df[(df['ref']=='C')]
            elif strand_assign == '-':
                df = df[(df['ref']=='G')]
                
        elif edit_type == 'ag':
            if strand_assign == '+':
                df = df[(df['ref']=='A')]
            elif strand_assign == '-':
                df = df[(df['ref']=='T')]
                
        pos_edits = pd.concat([pos_edits, df])

        l =[]

        for k,v in edits.items():
            for kv,vv in v.items():
                if vv>0:
                    l.append([chrom, int(kv), int(kv)+1, k, strand_assign])
        edits_df = pd.concat([edits_df, pd.DataFrame(l)])
    return pos_edits, edits_df

def assign_utr(x):
    strand = x[6]
    utr_start = x[3]
    utr_end = x[4]
    if strand =='+':
        if utr_end < gtf_CDS_start[x['isoform']]:
            return '5UTR'
        elif utr_start > gtf_CDS_end[x['isoform']]:
            return '3UTR'

    if strand =='-':
        if utr_start > gtf_CDS_end[x['isoform']]:
            return '5UTR'
        elif utr_end <  gtf_CDS_start[x['isoform']]:
            return '3UTR'
def parallelize_dataframe(df, func, n_cores=cpu_count()-1):
    '''Multiprocess getting edits for different regions within an isoform group'''
    df_split = np.array_split(df, n_cores)
    pool = Pool(n_cores)
    results = pool.map(func, df_split)
    pool.close()
    pool.join()
    
    # Concatenate the first object from each process
    n_positions = pd.concat([r[0] for r in results], axis=0)
    
    # Concatenate the second object from each process
    n_edited = pd.concat([r[1] for r in results], axis=0)
    
    return n_positions, n_edited

def filter_for_edits(x):
    if edit_type =='ct':
        if (x['strand']=='+') and (x['ref'] == 'C') and (x['T'] > 0):
            return x['T']
        elif (x['strand']=='-') and (x['ref'] == 'G') and (x['A'] > 0):
            return x['A']
        else:
            return 0
    elif edit_type =='ag':
        if (x['strand']=='+') and (x['ref'] == 'C') and (x['T'] > 0):
            return x['T']
        elif (x['strand']=='-') and (x['ref'] == 'G') and (x['A'] > 0):
            return x['A']
        else:
            return 0
        
def format_and_output_tables(table, sample):
    table['conversion'] = table.apply(filter_for_edits, axis=1)
    table['coverage'] = table.apply(lambda x: x['A'] + x['T'] + x['C'] + x['G'], axis=1)
    table = table.reset_index().rename(columns={'index' : 'start' } )
    table['end'] = table['start'] +1
    table['confidence'] = table.apply(lambda x: 1 - betainc(x['conversion'], (x['coverage'] - x['conversion']), 0.05), axis=1)
    table.to_csv(os.path.join(output_dir, sample + '.tsv'), sep='\t', index=False)
    
    table_edited = table[table['conversion']!=0]
    table_edited['key'] = table_edited.apply(lambda x: f"{x['gene']}:{x['isoform']}:{x['chrom']}:{x['start']}:{x['strand']}", axis=1)
    table_edited['edit_info'] = table_edited.apply(lambda x: f"{x['conversion']},{x['coverage']}", axis=1)
    table_edited = table_edited[['chrom', 'start', 'end', 'key', 'edit_info', 'strand']]
    table_edited.to_csv(os.path.join(output_dir, sample + '.bed'), sep='\t', index=False, header=None)



    
    
# Specifying arguments
parser = argparse.ArgumentParser(description='Identify read level edits. Make sure all split BAMs have index files prior to running script.')
parser.add_argument("--in_dir", type=str, default='./', help='directory containing all split BAM files', required=True)
parser.add_argument("--out_dir", type=str, default='./', help='directory where BED files will be saved', required=True)
parser.add_argument("--gtf_path", type=str, help='full path to GFF generated from HQ clustered PacBio reads and SQANTI3', required=True)
parser.add_argument("--ref", type=str, help='full path to reference fasta (e.g. hg19.fa)', required=True)
parser.add_argument("--split_pickle", type=str, help='full path to split bam dictionary pickle file (output from split bam script)', required=True)
parser.add_argument("--split", type=str, help='name of split (format: split_bam_#_#). these should match keys in the split pickle dictionary.', required=True)


args = parser.parse_args()
input_dir = args.in_dir
output_dir = args.out_dir
#read_level_out = args.read_edits_out_dir
#iso_counts_file = args.iso_counts
gtf = args.gtf_path
reffile = args.ref
split_pickle_path = args.split_pickle
split_name = args.split

sorted_bam = os.path.join(output_dir, split_name + '.sorted.bam')
#rbp_test_read_df = pd.read_csv(os.path.join(output_dir, iso_counts_file), sep='\t')

## Read in pickle file (split bam dictionary)
infile = open(split_pickle_path,'rb')
split_dictionary = pickle.load(infile)
infile.close()    

edit_type = 'ct'

# Make Transcript level GTF
print('Constructing transcript-level annotation information...')
gtf_transcript = pd.read_csv(gtf, sep='\t', comment='#', header=None, low_memory=False)
gtf_transcript = gtf_transcript[gtf_transcript[2] == 'transcript']
gtf_transcript['isoform'] = gtf_transcript[8].str.extract(r'transcript_id "([^"]+)"', expand=False)
gtf_transcript['gene'] = gtf_transcript[8].str.extract(r'gene_id "([^"]+)"', expand=False)

# Make Transcript CDS level GTF
gtf_CDS = pd.read_csv(gtf, sep='\t', comment='#', header=None)
gtf_CDS = gtf_CDS[gtf_CDS[2] == 'CDS']
gtf_CDS['isoform'] = gtf_CDS[8].str.extract('transcript_id "(.*?)";', expand=False)
gtf_CDS['gene'] = gtf_CDS[8].str.extract('gene_id "(.*?)";', expand=False)

gtf_CDS_start = dict(zip(gtf_CDS['isoform'], gtf_CDS[3]))
gtf_CDS_end = dict(zip(gtf_CDS['isoform'], gtf_CDS[4]))

# Make Transcript UTR level GTF and classify them
gtf_utr = pd.read_csv(gtf, sep='\t', comment='#', header=None)
gtf_utr = gtf_utr[gtf_utr[2] == 'UTR']
gtf_utr['isoform'] = gtf_utr[8].str.extract('transcript_id "(.*?)";', expand=False)
gtf_utr['gene'] = gtf_utr[8].str.extract('gene_id "(.*?)";', expand=False)
gtf_utr[2] = gtf_utr.apply(assign_utr, axis=1)
gtf_3utr = gtf_utr[gtf_utr[2]=='3UTR']
gtf_5utr = gtf_utr[gtf_utr[2]=='5UTR']

# Get GTF information only for the isoforms in the split
group_transcripts_transcript = gtf_transcript[gtf_transcript.isoform.isin(split_dictionary[split_name]['isoform'].tolist())]
group_transcripts_CDS = gtf_CDS[gtf_CDS.isoform.isin(split_dictionary[split_name]['isoform'].tolist())]
group_transcripts_5utr = gtf_5utr[gtf_5utr.isoform.isin(split_dictionary[split_name]['isoform'].tolist())]
group_transcripts_3utr = gtf_3utr[gtf_3utr.isoform.isin(split_dictionary[split_name]['isoform'].tolist())]

# Check that none of the entries have the same start and stop
group_transcripts_transcript = group_transcripts_transcript[group_transcripts_transcript[3]!= group_transcripts_transcript[4]]
group_transcripts_CDS = group_transcripts_CDS[group_transcripts_CDS[3]!= group_transcripts_CDS[4]]
group_transcripts_5utr = group_transcripts_5utr[group_transcripts_5utr[3]!= group_transcripts_5utr[4]]
group_transcripts_3utr = group_transcripts_3utr[group_transcripts_3utr[3]!=group_transcripts_3utr[4]]

# Get the position level and read level data for full transcript,CDS, 5', and 3' regions
print('Identifying edits...')
n_positions_transcript, n_edited_transcript = parallelize_dataframe(group_transcripts_transcript, get_position_matrix)
print('(1/4) Transcript')
n_positions_CDS, n_edited_CDS = parallelize_dataframe(group_transcripts_CDS, get_position_matrix)
print('(2/4) CDS')
n_positions_5utr, n_edited_5utr = parallelize_dataframe(group_transcripts_5utr, get_position_matrix)
print('(3/4) 5UTR')
n_positions_3utr, n_edited_3utr = parallelize_dataframe(group_transcripts_3utr, get_position_matrix)
print('(4/4) 3UTR')

read_level_out  = os.path.join(output_dir, f'read_level_edits_{edit_type}/')

# Prior to output-ing files, we check that the directory exists
if not os.path.exists(read_level_out):
    os.makedirs(read_level_out)
    print('making directory')

print('Creating output files...')

### CREATE FUNCTION TO MAKE THESE INTO A BED FILE
n_edited_transcript.to_csv(os.path.join(read_level_out, f"{split_name}_read_level_edits.bed"), header=None, sep='\t', index=False)
n_edited_CDS.to_csv(os.path.join(read_level_out, f"{split_name}_read_level_edits_CDS.bed"), header=None, sep='\t', index=False)
n_edited_5utr.to_csv(os.path.join(read_level_out, f"{split_name}_read_level_edits_5utr.bed"), header=None, sep='\t', index=False) 
n_edited_3utr.to_csv(os.path.join(read_level_out, f"{split_name}_read_level_edits_3utr.bed"), header=None, sep='\t', index=False)

if n_positions_transcript.shape[0] >0:
    format_and_output_tables(n_positions_transcript, f"{split_name}_{edit_type}_transcript")
if n_positions_CDS.shape[0] >0:
    format_and_output_tables(n_positions_CDS, f"{split_name}_{edit_type}_CDS")
if n_positions_5utr.shape[0] >0:
    format_and_output_tables(n_positions_5utr, f"{split_name}_{edit_type}_5utr")
if n_positions_3utr.shape[0] >0:
    format_and_output_tables(n_positions_3utr, f"{split_name}_{edit_type}_3utr")

pybedtools.helpers.cleanup(verbose=True,remove_all=True)

with open(os.path.join(output_dir, 'completed.log'), "a") as log:
    log.write(split_name + '\n')
#### in split bam file, create a log and append the split name to it