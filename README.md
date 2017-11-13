# Segmentum
Segmentum: a tool for copy number analysis of cancer genome.

## Introduction

Segmentum is a tool for the identification of CNAs and copy-neutral loss of heterozygosity (LOH) in tumor samples using whole-genome sequencing data. Segmentum segments the genome by analyzing the read-depth and B-allele fraction profiles using a double sliding window method. It requires a matched normal sample to correct for biases and to discriminate somatic from germline events. Segmentum, written in the Python programming language, is fast and performs segmentation of a whole genome in less than two minutes.

### Citation information
Scientific paper corresponding to this work is openly accessible online at BMC Bioinformatics.

Afyounian, E., Annala, M. and Nykter, M., 2017. Segmentum: a tool for copy number analysis of cancer genomes. BMC Bioinformatics, [online] 18(215), pp.1-10. Available at: <https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-017-1626-8>


## Installation
To install Segmentum, download the latest release and extract it to a folder. Then run the Makefile and add the subfolder bin/ to your PATH. See below for an example:    

```
wget https://github.com/eafyounian/Segmentum/archive/master.zip   
unzip master.zip   
cd Segmentum-master/    
make   
export PATH=/some/folder/Segmentum-master/bin:$PATH   
```
Segmentum requires [Samtools](https://github.com/samtools/samtools) to be installed.     

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
    
    CHROM   POSITION    REF     ALT     sample_x_t
    chr1    10174       C       T       0.00
    chr1    10177       A       C       0.08
    chr1    10583       G       A       0.00
    chr1    15211       T       G       0.71
    chr1    16487       T       C       0.15

## Output files
Segmentum outputs a SEG file containing the identified segments. 

**Example**  

    'ID              chrom   loc.start       loc.end     num_probes      seg.mean        baf.mean
    sample_x_t       chr1    8000            80000       72              0.246188        0.425905
    sample_x_t       chr1    80000           173000      93              -0.007791       0.365433
    sample_x_t       chr1    173000          270000      97              -0.708895       0.487773
     
## Running Segmentum
Segmentum can be run in two modes with/without using B-allele fraction data with the following commands:

```    
## with BAF
Segmentum analyze with BAF <tumor> <normal> <BAF_file> <window_size> <clogr_threshold> <BAF_threshold> [-m N] [-l N] [-b N] [-p N] [-B N]     

## without BAF
Segmentum analyze without BAF <tumor> <normal> <window_size> <clogr_threshold> [-m N] [-l N]  


Options:  
    -h --help             Show this screen
    -m --min_read=N       Minimum number of reads from the normal sample to calculate the coverage log ratio [default: 50]  
    -l --logr_merge=N     Log ratio segments merging threshold [default: 0.15]  
    -b --baf_merge=N      B-allele fraction segments merging threshold [default: 0.05] 
    -p --print=N          If true, prints the results to standard output, otherwise to a file with the same name as the sample name [default: True]    
    -B --BAF=N            If true, creates a .WIG file for heterozygous SNP allele fractions to be opened in IGV [default: True]    
    
```

**Example with BAF**   
    `Segmentum analyze with BAF sample_x_t.wig.gz sample_x_n.wig.gz B_allele_fraction.tsv.gz 11 0.7 0.3`

Segmentum outputs the results to standard output by default. Use `-p False`  to create a SEG file with the same name as the tumor file.
    
## Visualization of the results in IGV
In order to visualize the results in IGV, a new file should be created containing only first 6 fields of the output file.

**Example**  
    `cut -f1,2,3,4,5,6 sample_x_t.seg  > IGV_sample_x_t.seg`  

`IGV_sample_x_t.seg` is now ready to be loaded in IGV.
    
## Extracting copy neutral LOH regions
In order to extract copy neutral LOH regions from the output(s) use the following command:

```
Segmentum find recurrent cnLOHs <seg_files>... [-c N] [-t N]  
    
Options:  
    -c --clogr_thresh=N   Coverage logratio must be below this threshold to call a copy neutral LOH region [default: 0.1]    
    -t --baf_thresh=N     B-allele fraction must be below this threshold to call a copy neutral LOH region [default: 0.15]   
```

**Example**  
`Segmentum find recurrent cnLOHs sample_x_t.seg`  
`Segmentum find recurrent cnLOHs sample_x_t.seg sample_y_t.seg sample_z_t.seg  #in case of more samples`      

## Creating the input files from BAM files
In order to create the coverage files, use the following command:

```    
Segmentum extract read depth <BAM_file> <window_size> [-q N]    
        
Options:   
    -q --quality=N        Minimum mapping quality [default: 10]   
```
    
**Example**  
```
Segmentum extract read depth sample_x_t.bam 2000 | gzip -c > sample_x_t.wig.gz
Segmentum extract read depth sample_x_n.bam 2000 | gzip -c > sample_x_n.wig.gz
```
        
In order to create the B-allele-fraction file, use the following command:

```    
Segmentum calculate BAF <genome_fasta> <SNP_position> <tumor> <normal> [--hetz=N:R] [-q N] [-r REGION]   
        
Options:  
    -r <region>           Restrict analysis to chromosomal region     
    -q --quality=N        Minimum mapping quality [default: 10]     
    --hetz=N:R            Minimum evidence for heterozygous [default: 4:0.3]    
```

**Example**  
```
Segmentum calculate BAF hg19.fa hg19_1000g_2014oct_SNPs.tsv.gz sample_x_t.bam sample_x_n.bam --hetz=4:0.3 -q20 | gzip -c > B_allele_fraction.tsv.gz
```

**Note on SNP postion file:** You can either use the provided SNP postion file or make one of your own. If you are planning to use the provided one you can download it from here: http://compbio.uta.fi/segmentum/hg19_1000g_2014oct_SNPs.tsv.gz. If you are planning to make your own SNP position file, the file should have the following format: 

```
    chr     position
    1       10177  
    1       10235  
    1       10352  
    1       10352  
    1       10505  

```

**Samples** You can find a set of sample coverage and BAF files for experimentation with Segmenum here:  http://compbio.uta.fi/segmentum/sample_files.tar.gz.    
     
**simulated data** You also can find a set of simulated coverage and BAF data as well as its corresponding ground truth for experimentation here: http://compbio.uta.fi/segmentum/simulated_data_example.tar.gz.    
