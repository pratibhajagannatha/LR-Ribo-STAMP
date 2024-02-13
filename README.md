# LR-Ribo-STAMP
Repository containing scripts to process long-read (LR) Ribo-STAMP data. 


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

```isoquant.py --reference /path/to/LR-Ribo-STAMP/test/ref/hg19_chr1.fa --genedb /path/to/LR-Ribo-STAMP/test/ref/gencode.v19.annotation.chr1.gtf --bam /path/to/LR-Ribo-STAMP/test/input_dir/Ribo-STAMP_HEK293T_rep1_chr1.readfiltered.sorted.bam --data_type pacbio --prefix Ribo-STAMP_HEK293T_rep1_chr1 --transcript_quantification unique_only --gene_quantification unique_only --output /path/to/LR-Ribo-STAMP/test/output_dir/ --threads 32```






















