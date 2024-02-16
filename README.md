# LR-Ribo-STAMP
Repository containing scripts to process long-read (LR) Ribo-STAMP data. 

## About LR-Ribo-STAMP
LR-Ribo-STAMP is an integrated experimental and computational method that couples long read sequencing with Ribo-STAMP to simultaneously measure transcription and translation at gene and mRNA isoform resolution. 

![main_schematic](https://github.com/pratibhajagannatha/LR-Ribo-STAMP/blob/main/LR-ribostamp_github_main.png)

## Suggested folder structure

```
├── sample_name
│   ├── input_dir
│   │   ├── aligned long reads (BAM)
│   │   ├── filtered aligned long reads (BAM)
│   ├── output_dir
│   │   ├── qc
│   │   ├── isoquant_out
│   │   ├── read_level_edits_ct
│   │   ├── split bam files and indices (.bam, .bam.bai)
│   │   ├── split bam dictionary
│   ├── ref
│   │   ├── fasta file
│   │   ├── GTF file
```

## Running the Computational Pipeline

### Step 1
QC the aligned long reads using [NanoPlot](https://github.com/wdecoster/NanoPlot). Although ```--tsv_stats``` and ```--raw``` are not required for running the tool, the outputs will be used for filtering aligned reads in Step 2.

```NanoPlot --raw --tsv_stats -o /path/to/LR-Ribo-STAMP/test/output_dir/qc/Ribo-STAMP_HEK293T_rep1_chr1/ --bam /path/to/LR-Ribo-STAMP/test/input_dir/Ribo-STAMP_HEK293T_rep1_chr1.sorted.bam```

### Step 2
Filter the aligned reads to remove low quality, non-uniquely mapped, and supplementary alignment reads. Reads mapped to the incorrect strand are also removed. After filtering, the script will outline the number of reads filtered out in total and by criteria. 

```python3 /path/to/LR-Ribo-STAMP/scripts/filter_bam_v2.py --out_dir /path/to/LR-Ribo-STAMP/test/input_dir/ --bam_path /path/to/LR-Ribo-STAMP/test/input_dir/Ribo-STAMP_HEK293T_rep1_chr1.sorted.bam --gtf /path/to/LR-Ribo-STAMP/test/ref/gencode.v19.annotation.chr1.gtf --qc_dir /path/to/LR-Ribo-STAMP/test/output_dir/qc/Ribo-STAMP_HEK293T_rep1_chr1/ --qc_cutoff 20```

Example output:
```
Number of reads filtered based on Q score: 16
Number of reads filtered based on mismatched strand: 568
Number of reads filtered based on supplementary read flag:578
Total reads filtered out: 1106
```

### Step 3
Run [IsoQuant](https://github.com/ablab/IsoQuant) on the filtered, aligned reads to obtain gene and isoform counts along with read assignments. 

```isoquant.py --reference /path/to/LR-Ribo-STAMP/test/ref/hg19_chr1.fa --genedb /path/to/LR-Ribo-STAMP/test/ref/gencode.v19.annotation.chr1.gtf --bam /path/to/LR-Ribo-STAMP/test/input_dir/Ribo-STAMP_HEK293T_rep1_chr1.readfiltered.sorted.bam --data_type pacbio --prefix Ribo-STAMP_HEK293T_rep1_chr1 --transcript_quantification unique_only --gene_quantification unique_only --output /path/to/LR-Ribo-STAMP/test/output_dir/Ribo-STAMP_HEK293T_rep1_chr1_isoquant --threads 32```

### Step 4
Split the bam file into smaller bam file for easy mutliprocessing and calling edits at mRNA isoform level. To create the split BAM files, each isoform of a gene having reads is assigned a random number (e.g. Isoform 1, Isoform 2, ..., Isoform N). The reads assigned to isoforms are then placed into the appropriate split BAM based on the isoforms' assigned number. This step also outputs a split_bam dictionary (in [pickle](https://docs.python.org/3/library/pickle.html) file format) that can be used to map isoforms and reads back to splits.

![split_bam_diagram](https://github.com/pratibhajagannatha/LR-Ribo-STAMP/blob/main/LR-ribostamp_splitbam_github_fig.png)

```python3 /path/to/LR-Ribo-STAMP/scripts/split_bam_isoquant.py --isoquant_dir /path/to/LR-Ribo-STAMP/test/output_dir/Ribo-STAMP_HEK293T_rep1_chr1_isoquant/ --genome_bam /path/to/LR-Ribo-STAMP/test/input_dir/Ribo-STAMP_HEK293T_rep1_chr1.readfiltered.sorted.bam --output_dir /path/to/LR-Ribo-STAMP/test/output_dir/ --sample_name Ribo-STAMP_HEK293T_rep1_chr1```


### Step 5
Use the make_read_level_script_ct.py script to create a bash script that will run edit calling on all of the split BAM files created in the previous step. If you are running this analysis on a server with a job scheduler/resource management system, you can run this command with ```--submitter``` to create a bash script that will submit an array job. The header may need to be modified for your scheduler (e.g. SLURM, TORQUE). Using this option will speed up the runtime significantly. 

```python3 /path/to/LR-Ribo-STAMP/scripts/make_read_level_edit_script_ct.py --working_directory /path/to/LR-Ribo-STAMP/test/output_dir/ --gtf_path /path/to/LR-Ribo-STAMP/test/ref/gencode.v19.annotation.chr1.gtf --split_pickle /path/to/LR-Ribo-STAMP/test/output_dir/split_bam_Ribo-STAMP_HEK293T_rep1_chr1 --ref /path/to/LR-Ribo-STAMP/test/ref/hg19.fa --job Ribo-STAMP_HEK293T_rep1_chr1 --script_path /path/to/LR-Ribo-STAMP/scripts/read_level_quant_se_ct_annotated.py```


### Step 6

#### A. Creating Gene and Isoform Coordinate FASTA files

If you do not already have gene and isoform coordinate fasta files, create a BED file from your GTF for only the exons and UTRs for genes or isoforms. There are a variety of tools you can use to do this. Once you have the BED file, sort the BED file using the following:

```sort -k1,1 -k2,2n output_exons_utrs.transcript.bed > output_exons_utrs.transcript.sorted.bed```

Merge overlapping exons and UTRs using the following:

```bedtools merge -i output_exons_utrs.transcript.sorted.bed -s -c 4,6 -o distinct > merged_exons_utrs.transcript.hg19.bed```

Obtain the FASTA sequences for the resulting coordinates:

```bedtools getfasta -fi hg19.fa -bed merged_exons_utrs.transcript.hg19.bed -name -s > merged_exons_utrs.transcript.hg19.fasta```

#### B. Filter Edits and Calculate EditsC

Once edits have been called successfully on all the split BAM files, the filter_edits_calc_editsC.py script can be used to filter edits and calculate the EditsC metric at gene and isoform levels. The script requires a sample file that contains the sample name, replicate number, condition, and path to the parent directory (sample_name in the folder structure outlined above). The script assumes that the necessary files are within the ```output_dir``` directory. The sample file should look like the following:

```
sample_name     rep     condition       path
APO1_1  rep1    APO     /path/to/APO1_1
APO1_2  rep2    APO     /path/to/APO1_2
APO1_3  rep3    APO     /path/to/APO1_3
RPS2_1  rep1    RPS2    /path/to/RPS2_1
RPS2_2  rep2    RPS2    /path/to/RPS2_2
RPS2_3  rep3    RPS2    /path/to/RPS2_3
```


Once the sample file is created, run the script similar to the following:

```python3 /path/to/LR-Ribo-STAMP/scripts/filter_edits_calc_editsC.py --sample_file /path/to/sample_file.txt --snp_file /path/to/hg19.commonSNPs147.bed3 --gtf_file/path/to/gencode.v19.annotation.gtf --gene_coordinates_fasta /path/to/merged_exons_utrs.genes.hg19.fasta --transcript_coordinates_fasta /path/to/merged_exons_utrs.transcript.hg19.fasta --output_dir /path/to/output```




















