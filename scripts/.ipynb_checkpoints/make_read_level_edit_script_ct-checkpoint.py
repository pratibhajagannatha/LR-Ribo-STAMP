#!/usr/bin/env python

import os
import pandas as pd
import pickle
from qtools import Submitter
import argparse

"""
Author: Pratibha Jagannatha (pjaganna@ucsd.edu)
Input files: working directory (generally where split BAM files are located), GTF, path to split dictionary (pickle file), reference Fasta, job name, path to read_level_quant_se_ct_annotated.py
Output file: bash script

Script to make the script to call edits on all split BAM files. 

Example: python3 /path/to/make_read_level_edit_script_ct.py --working_directory /path/to/output_dir/ --gtf_path /path/to/refs/gencode.v19.annotation.gtf --split_pickle /path/to/split_bam_dictionary --ref /path/to/refs/hg19.fa --job job_name --script_path /path/to/read_level_quant_se_ct_annotated.py

"""

parser = argparse.ArgumentParser(description='Split the bam file into chunks that can be processed in parallel.')
#parser.add_argument("--working_dir", type=str, default='./', help='where cluster report file and group.txt is located', required=True)
parser.add_argument("--working_directory", type=str, default='./', help='full path to directory where split bams are and output directory will be made', required=True)
parser.add_argument("--gtf_path", type=str, help='full path to PacBio GFF file (unannotated) or Gencode GTF annotation (annotated)', required=True)
parser.add_argument("--split_pickle", type=str, help='full path to split bam dictionary pickle file', required=True)
parser.add_argument("--ref", type=str, help='full path to reference annotation', required=True)
parser.add_argument("--job", type=str, help='name for jub submission', required=True)
parser.add_argument("--script_path", type=str, default='./read_level_quant_se_ct.py', help='full path to read_level_quant_se_ct_annotated.py script', required=True)
parser.add_argument("--submitter", help="increase output verbosity", action="store_true")
args = parser.parse_args()
working_dir = args.working_directory
cluster = args.gtf_path
pickle_path = args.split_pickle
ref_path = args.ref
job_name = args.job
script = args.script_path

infile = open(pickle_path,'rb')
split_bam_dict = pickle.load(infile)
infile.close()

isExist = os.path.exists('{}/read_level_edits_ct'.format(working_dir))
if not isExist:
    os.makedirs('{}/read_level_edits_ct'.format(working_dir))
    print("Output directory does not exist, making new directory.")

print('There are {} items in the dictionary.'.format(len(split_bam_dict)) + '\n')

print('Constructing commands...\n')

cmds = []
for k in list(split_bam_dict.keys()):
    cmd = 'python3 {}'.format(script) 
    cmd += ' --in_dir {}'.format(working_dir)
    cmd += ' --out_dir {}'.format(working_dir)
    cmd += ' --gtf_path {}'.format(cluster)
    cmd += ' --ref {}'.format(ref_path)
    cmd += ' --split_pickle {}'.format(pickle_path)
    cmd += ' --split {}'.format(k)
    cmds.append(cmd)

print('Exaple command : '  + cmds[1] + '\n')

print('Making script...\n')

if args.submitter:
    Submitter(
        commands=cmds,
        job_name=job_name,
        sh=os.path.join(working_dir,'read_level_edits_ct/{}.sh'.format(job_name)),
        array=True,
        nodes=1,
        ppn=4,
        walltime='12:00:00',
        submit=False,
    )
else:
    with open(os.path.join(working_dir,'read_level_edits_ct/{}.sh'.format(job_name)), 'w') as cmd_writer:
        cmd_writer.write('#!/bin/bash\n')
        for c in cmds:
            cmd_writer.write(c + '\n')