# Segmentum
Segmentum a fast tool for copy number analysis of cancer genome.

##Introduction

Segmentum is a fast tool for the identification of CNAs and copy-neutral loss of heterozygosity (LOH) in tumor samples using whole-genome sequencing data. Segmentum segments the genome by analyzing the read-depth and B-allele fraction profiles using a double sliding window method. It requires a matched normal sample to correct for biases such as GC-content and mapability and to discriminate somatic from germline events. Segmentum, written in the Python programming language, is fast and performs segmentation of a whole genome in less than two minutes.
    
## Input files
**1.** Tumor and its matched normal coverage files with WIG format.
        
**Example**  
        
        fixedStep chrom=chr20 start=1 step=1000
        140
        151
        159
        144
        131

**2.** B-allele-fraction file containing the B-allele fractions at heterozygous SNPs.

**Example**  
    
    CHROM   POSITION    REF     ALT     G9_6338_t
    chr1    10174       C       T       0.00
    chr1    10177       A       C       0.08
    chr1    10583       G       A       0.00
    chr1    15211       T       G       0.71
    chr1    16487       T       C       0.15

## Output files
Segmentum outputs a SEG file containing the identified segments. 
    
**Example**  

    'ID             chrom   loc.start       loc.end         seg.mean        baf.mean
    G9_6338_t       chr14   22425001        23000001        0.246188        0.425905
    G9_6338_t       chr14   23000001        86937001        -0.007791       0.365433
    G9_6338_t       chr14   86937001        88319001        -0.708895       0.487773
    G9_6338_t       chr14   88319001        88473001        0.016230        0.433023
    G9_6338_t       chr14   88473001        88926001        -0.672465       0.456555
    G9_6338_t       chr14   88926001        107289001       0.022306        0.328848
     
## Running Segmentum
Segmentum can be run with the following command:

```    
Segmentum.py <tumor_coverage> <normal_coverage> <B_allele_fraction> <window_size> <clogr_threshold> <BAF_threshold> [-l N] [-b N] [-m N] [-z N]  

Options:  
    -h --help             Show this screen  
    -l --logr_merge=N     Logratio merging threshold [default: 0.15]  
    -b --baf_merge=N      B-allele fraction merging threshold [default: 0.05]  
    -m --min_read=N       Minimum number of reads for a position to be considered while calculating the coverage logratio [default: 50]  
    -z --z_score=N        Number of standard deviations away from the mean to call a breakpoint [default: 5]    
```

**Example**   
    `python Segmentum.py G9_6338_t.wig.gz G9_6338_n.wig.gz B_allele_fraction.tsv.gz 11 0.8 0.3`

Segmentum outputs a SEG file with the same name as the tumor file.
    
## Visualization of the results in IGV
In order to visualize the results in IGV a new file should be created containing only first 5 fields of the output file.

**Example**  
    `cut -f1,2,3,4,5 G9_6338_t.seg  > IGV_g9_6338.seg`
The IGV_g9_6338.seg is now ready to be loaded in IGV.
    
## Extracting copy neutral LOH regions
In order to extract copy neutral LOH regions from the output(s) following command is used:

```
python Recurrent_cnLOH.py <seg_files>... [-c N] [-b N]  
    
Options:  
    -h --help         Show this screen.
    -c --clogr_thresh=N   Coverage logratio must be below this threshold to call a copy neutral LOH region [default: 0.1]
    -b --baf_thresh=N     B-allele fraction must be below this threshold to call a copy neutral LOH region [default: 0.15]  
```

**Example**  
    
    `python Recurrent_cnLOH.py G9_6338_t.seg`

## Creating the input files from BAM files
In order to create the coverage files, one can use the Pypette package available at: (https://github.com/annalam/pypette) using the following command:

```    
coverage tiled <bam_file> <window_size> [-s N] [-q N] [-S|-1|-2] [-P|-M]  
        
Options:   
    -q --quality=N Minimum alignment quality [default: 10].  
    -s --step=N Step size for window placement [default: window size / 2].   
    -S --single Use all reads for coverage calculation, not just paired.   
    -P --plus Calculate coverage only for the plus strand.   
    -M --minus Calculate coverage only for the minus strand.  
```
    
**Example**  
    
    ```
    coverage tiled G9_6338_t.bam 2000 | gzip -c > G9_6338_t.wig.gz
    coverage tiled G9_6338_n.bam 2000 | gzip -c > G9_6338_n.wig.gz
    ```
        
In order to create the B-allele-fraction file, one can use the Pypette package available at: (https://github.com/annalam/pypette) using the following commands:

```    
variant call <genome_fasta> <bam_files>... [-r REGION] [--ref=N:R] [--hetz=N:R] [--homz=N:R] [-q N] [-Q SAMPLES] [--keep-all]  
variant keep samples <vcf_file> <regex>  
variant heterozygous bases <vcf_file> <pos_file>   
variant discard samples <vcf_file> <regex>  
variant allele fractions <vcf_file> <pos_file>   
        
Options:  
    -r <region> Restrict analysis to chromosomal region  
    -q N Minimum mapping quality score [default: 10]  
    -Q SAMPLES Samples for which mapping quality is ignored [default: ]  
    --ref=N:R Minimum evidence for homozygous reference [default: 8:0.9]  
    --hetz=N:R Minimum evidence for heterozygous [default: 4:0.25]  
    --homz=N:R Minimum evidence for homozygous alt [default: 4:0.8]  
    --keep-all Show sites even if they are all homozygous reference    
```

**Example**  
        
1. Run the following command for each chromosome (example command for chr20 is shown):
    `variant call --hetz=4:0.15 --homz=4:0.8 --ref=8:0.95 -Q _n -r chr20 hg19.fa G9_6338_t.bam G9_6338_n.bam | gzip -c > chr20.vcf.gz`
2. `cat <(gunzip -c chr1.vcf.gz | grep 'CHROM') <(cat chr{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y}.vcf.gz | gunzip -c | grep -v 'CHROM') | pigz -c > variants.vcf.gz`
3. `variant keep samples variants.vcf.gz G9_6338_n.bam | variant heterozygous bases - heterozygous_SNPs_position_file > heterozygous_snps.tsv`
4. `gzip heterozygous_snps.tsv`
5. `variant discard samples variants.vcf.gz G9_6338_n.bam | variant allele fractions - heterozygous.tsv.gz | gzip -c > B_allele_fractions.tsv.gz`
