#!/bin/env python

"""
A tool for copy number analysis and segmenting the cancer genome.

Usage:
  segmentum extract read depth <BAM_file> <window_size> [-q N] 
  segmentum calculate BAF <genome_fasta> <SNP_position> <tumor> <normal> [--hetz=N:R] [-q N] [-r REGION]
  segmentum analyze with BAF <tumor> <normal> <BAF_file> <window_size> <clogr_threshold> <BAF_threshold> [-m N] [-l N] [-b N] [-p N] [-B N]
  segmentum analyze without BAF <tumor> <normal> <window_size> <clogr_threshold> [-m N] [-l N] [-p N]
  segmentum find recurrent cnLOHs <seg_files>... [-c N] [-t N]
  segmentum simulate <normal> [-P N] [-O N] [-L N]
  
Options:
  -h --help             Show this screen.
  -m --min_read=N       Minimum number of reads from the normal sample to calculate the coverage log ratio [default: 50]
  -l --logr_merge=N     Log ratio segments merging threshold [default: 0.15]
  -b --baf_merge=N      B-allele fraction segments merging threshold [default: 0.05]
  -q --quality=N        Minimum mapping quality [default: 10]
  -r <region>           Restrict analysis to chromosomal region
  --hetz=N:R            Minimum evidence for heterozygous [default: 4:0.3]
  -c --clogr_thresh=N   Coverage logratio must be below this threshold to call a copy neutral LOH region [default: 0.1]
  -t --baf_thresh=N     B-allele fraction must be below this threshold to call a copy neutral LOH region [default: 0.15]
  -p --print=N          If true, prints the results to standard output, otherwise to a file with the same name as the sample name [default: True]
  -B --BAF=N            If true, creates a .WIG file for heterozygous SNP allele fractions to be opened in IGV [default: True]
  -P --tumor_purity=N   (1 - fraction) of normal contamination [default: 0.7]
  -O --output_prefix=N  prefix to be assigned to the simulated files [default: simulated]
  -L --read_length=N    read length to be considered for simulation [default: 150]

Author: Ebrahim Afyounian <ebrahim.afyounian@uta.fi>
"""

from __future__ import print_function
from coverage_tiled import coverage_tiled
from variant_call import calculate_BAF
from Recurrent_cnLOH import find_chromosome_length, find_cnloh, find_recurrent
from simulator import simulate, open_file
import sys, gzip, math, time, re, docopt, subprocess, os
import scipy.signal as sig 
import numpy as np

class Chromosome:
    """class to keep track of start, step, and (either coverage OR logRatios) data from a wig file"""
    def __init__(self):
        self.start = None 
        self.step = None
        self.values = []

class Options(object):
    def __init__(self, args):
        self.region = args['-r']
        self.min_mapq = int(args['--quality'])
        self.ignore_mapq = None
        self.keep_all = None
        self.min_ref_reads = 8
        self.min_ref_ratio = 0.9
        hetz = args['--hetz'].split(':')
        self.min_hetz_reads = int(hetz[0])
        self.min_hetz_ratio = float(hetz[1])
        self.min_homz_reads = 4
        self.min_homz_ratio = 0.7

class AlleleFraction:
    """class to keep track of positions and values of heterozygous SNPs"""
    def __init__(self):
        self.values = [] 
        self.positions = []

def read_wig(sample_file):
    """reads a wig file and returns a dictionary out of it"""
    if sample_file.endswith('.gz'):
        file = subprocess.Popen('gunzip -c %s' % sample_file, stdout=subprocess.PIPE, shell=True).stdout
    else:
        file = open(sample_file, 'r')

    chromosomes = {}
    chr = None
    for line in file:
        line = line.decode('utf8')
        if line[0] == 't': continue    ## checking whether line starts with 'track'
        if line[0] in 'fv':            ## checking whether line starts with 'fixedStep' or 'variableStep'
            chrom = re.search('chrom=(\S+)', line).group(1)
            start = int(re.search('start=(\d+)', line).group(1))
            step = int(re.search('step=(\d+)', line).group(1))
            
            chr = Chromosome()
            chr.start = start
            chr.step = step
            values = []
            chr.values = values
            
            chromosomes[chrom] = chr
        else:
            values.append(float(line))
    file.close()
    return chromosomes

def make_wig_file(allele_fractions, out_path):
    """makes a wig file out of mean_filtered_allele_fractions"""
    file = open(out_path,'w')
    file.write("track graphType=bar windowingFunction=none autoScale=off viewLimits=-2:2\n")
    for chrom in allele_fractions.keys():
        start = allele_fractions[chrom].start
        step = allele_fractions[chrom].step
        values = allele_fractions[chrom].values
        file.write("fixedStep chrom=%s start=%d step=%d\n" %(chrom, start, step))

        for value in values:
            if str(value) == 'nan':
                file.write('NaN' + '\n')
            else:
                file.write(str(value) + '\n')
    file.close()

def make_wig_file2(dict, outPath):
    """makes a wig file out of allele fraction data. Receives a dictionary"""
    file = open(outPath, 'w')
    file.write('track graphType=bar windowingFunction=none autoScale=off viewLimits=0:0.5\n')
    for chrom in dict.keys():
        values  = dict[chrom].values
        positions = dict[chrom].positions
        file.write("variableStep chrom=" + chrom + "\n")
        for val in zip(positions, values):
            file.write('%s %s\n' %(val[0], val[1] if str(val[1]) != 'nan' else 'NaN'))
    file.close()

def extract_sample(tsv_path, sample_name):
    """Extracts the allele fraction data based on sample name; Returns a dictionary"""
    if tsv_path.endswith('.gz'):
        file = subprocess.Popen('gunzip -c %s' % tsv_path, stdout=subprocess.PIPE, shell=True).stdout
    else:
        file = open(tsv_path, 'r')
    headers = next(file).decode('utf8').strip().split('\t')
    if 'gz' in sample_name: sample_name = sample_name.replace('.gz', '')
    if 'wig' in sample_name: sample_name = sample_name.replace('.wig', '')
    for idx, col in enumerate(headers):
        if re.search(sample_name, col) == None: continue
        else:
            sampleIdx = idx
            break

    allele_fractions = {}
    chrom = None

    for line in file:
        line = line.decode('utf8')
        if not line: continue 
        lst = line.split('\t')
        if lst[0] != chrom:
            chrom = lst[0]
            al_frac = AlleleFraction()
            
            values = []
            al_frac.values = values
            
            positions = []
            al_frac.positions = positions
            values.append( float(lst[sampleIdx]) )
            positions.append( int(lst[1]) )

            allele_fractions[chrom] = al_frac
        else:
            values.append( float(lst[sampleIdx]) )
            positions.append( int(lst[1]) )
    return allele_fractions

def smooth_allele_fractions(allele_fractions, win_size):
    """Smooths the allele fraction data and returns it"""
    for chrom in allele_fractions:
        allele_fracs = allele_fractions[chrom].values
        allele_fracs = np.array(allele_fracs) #####this added on 15th FEB 16
        # l1 = sig.medfilt(abs(np.subtract(0.5, allele_fracs)), win_size)  #deltaY for 0/100 allele fraction
        l1 = sig.medfilt(np.absolute(np.subtract(0.5, allele_fracs)), win_size)  #deltaY for 0/100 allele fraction
        # l2 = abs(0.5 - sig.medfilt(allele_fracs, win_size)) #deltaY for 50/50 allele fraction
        l2 = np.absolute(0.5 - sig.medfilt(allele_fracs, win_size)) #deltaY for 50/50 allele fraction
        # nearness = 1 - (2 * abs(np.subtract(0.5, allele_fracs)))
        nearness = 1 - (2 * np.absolute(np.subtract(0.5, allele_fracs)))
        allele_fractions[chrom].values = nearness * l2 + ((1 - nearness) * l1)
    return allele_fractions

def allele_frac_mean_filt(allele_fractions, win_size, step, start):
    """applies a mean filter over the allele fractions"""

    mean_filtered_allele_fractions = {}
    for chrom in allele_fractions:
        positions = allele_fractions[chrom].positions
        values = allele_fractions[chrom].values
        avg_allele_frac = []
        chr_size = positions[-1]
        win_start = start
        win_end = win_size * step
        first_snp = 0     ##index of the snp
        first_out_snp = 0
        win_sum = 0
        win_count = 0

        while win_start < chr_size:
            while positions[first_snp] < win_start:
                win_sum -= values[first_snp]
                win_count -= 1
                first_snp += 1
            while first_out_snp < len(positions) and positions[first_out_snp] <= win_end:
                win_sum += values[first_out_snp]
                win_count += 1
                first_out_snp += 1
            if win_count < (win_size/2) or win_count < 10:    ##if there is not as much snp as win_size, put nan #can be made as an option chosen by the user
                avg_allele_frac.append(np.nan) #nan
            else:        
                avg_allele_frac.append(win_sum / win_count)
            win_start += step
            win_end += step
        
        avg_allele_frac = np.concatenate((   np.full( int(win_size / 2)+1 , np.nan   )  , avg_allele_frac))
        
        chr = Chromosome()
        chr.start = start
        chr.step = step
        chr.values = avg_allele_frac
        mean_filtered_allele_fractions[chrom] = chr
    return mean_filtered_allele_fractions

def find_mean_difference(allele_fractions, win_size):
    """finds the absolute mean difference of the allele fractions for all the chromosomes"""
    mean_diffs = {}
    
    for chrom in allele_fractions:
        means = np.array(allele_fractions[chrom].values)
        abs_mean_diffs = np.absolute( means[(win_size - 1):] - means[: means.size - win_size + 1])
        abs_mean_diffs = np.concatenate((np.zeros(int((win_size - 1) / 2)), abs_mean_diffs ))
        abs_mean_diffs = np.concatenate((abs_mean_diffs , np.zeros(int((win_size - 1) / 2))))
        
        chr = Chromosome()
        chr.start = allele_fractions[chrom].start
        chr.step = allele_fractions[chrom].step
        chr.values = abs_mean_diffs
        mean_diffs[chrom] = chr
    return mean_diffs

def estimate_threshold(win_size, logr_thresh, baf_thresh, logr_merging_thresh, baf_merging_thresh, z_score):
    
    win_sizes = []
    logr_thresholds = []
    baf_thresholds = []
    
    final_win_size = int(round( win_size / ( float(logr_merging_thresh) / logr_thresh ) ** 2 ))
    if final_win_size % 2 == 0: final_win_size += 1
    
    ## Finding the logr thresholds
    while logr_thresh > logr_merging_thresh:
        win_sizes.append(win_size)
        logr_thresholds.append(round(logr_thresh, 2))
        standard_deviation = logr_thresh / z_score
        win_size = int(round(1.5 * win_size))
        if win_size % 2 == 0: win_size += 1
        logr_thresh = z_score * standard_deviation * (math.sqrt(2.0 / 3))
    
    win_sizes.append(final_win_size)
    logr_thresholds.append(logr_merging_thresh)
    
    ## Finding the baf thresholds
    baf_thresholds.append(baf_thresh)
    for win_size in range(1, len(win_sizes)):
        baf_thresh = round(baf_thresh * math.sqrt(2.0 / 3), 2)
        baf_thresholds.append(baf_thresh)
    
    for index, baf_thresh in enumerate(baf_thresholds):
        if baf_thresh < baf_merging_thresh:
            baf_thresholds[index] = baf_merging_thresh
    
    sys.stderr.write("The following window_sizes were estimated: ") 
    for size in win_sizes: sys.stderr.write(str(size) + ", ")
    sys.stderr.write("\nThe following thresholds were estimated for log ratios track: ") 
    for threshold in logr_thresholds: sys.stderr.write(str(threshold) + ", ")
    sys.stderr.write("\nThe following thresholds were estimated for BAF track: ") 
    for threshold in baf_thresholds: sys.stderr.write(str(threshold) + ", ")
    sys.stderr.write("\n") 
    return win_sizes, logr_thresholds, baf_thresholds

def estimate_threshold_only_for_coverage(win_size, logr_thresh, logr_merging_thresh, z_score):
    
    win_sizes = []
    logr_thresholds = []

    final_win_size = int(round( win_size / ( float(logr_merging_thresh) / logr_thresh ) ** 2 ))
    if final_win_size % 2 == 0: final_win_size += 1
    
    while logr_thresh > logr_merging_thresh:
        win_sizes.append(win_size)
        logr_thresholds.append(round(logr_thresh, 2))
        standard_deviation = logr_thresh / z_score
        win_size = int(round(1.5 * win_size))
        if win_size % 2 == 0: win_size += 1
        logr_thresh = z_score * standard_deviation * (math.sqrt(2.0 / 3))
        
        win_sizes.append(final_win_size)
        logr_thresholds.append(logr_merging_thresh)

    sys.stderr.write("The following window_sizes were estimated: ")
    for size in win_sizes: sys.stderr.write(str(size) + ", ")
    sys.stderr.write("\nThe following thresholds were estimated for log ratios track: ")
    for threshold in logr_thresholds: sys.stderr.write(str(threshold) + ", ")
    sys.stderr.write("\n")
    return win_sizes, logr_thresholds

def double_sliding_window_all_chromosomes(log_ratios, win_size):
    """finds the mean differences for the coverage log ratios of all the chromosomes"""
    
    mean_differences = {}
    for chrom in log_ratios:
        ## Using np.convolve() as a mean filter (Implements the mean filter using convolution)
        mean_filtered_values = np.convolve(log_ratios[chrom].values, np.ones((win_size,))/win_size, mode='same')  
        chr_mean_diffs = np.absolute(mean_filtered_values[(win_size - 1):] - mean_filtered_values[: mean_filtered_values.size - win_size + 1])
        chr_mean_diffs = np.concatenate((np.zeros(int((win_size - 1) / 2)), chr_mean_diffs))
        chr_mean_diffs = np.concatenate((chr_mean_diffs, np.zeros(int((win_size - 1) / 2))))
        
        chr = Chromosome()
        chr.start = log_ratios[chrom].start
        chr.step = log_ratios[chrom].step
        chr.values = chr_mean_diffs
        mean_differences[chrom] = chr
    return mean_differences

def filter_outliers(log_ratios, win_size):
    """filters out outliers in logRatios using the median filter and returns the log_ratios"""
    for chrom in log_ratios:
        log_ratios[chrom].values = sig.medfilt(log_ratios[chrom].values, win_size)
    return log_ratios

def calculate_logratios(tumor_sample, reference_sample, min_read):
    """Builds the dictionary that holds the logRatio for each chromosome and also takes care of the systematic bias"""
    log_ratios = {} #dictionary holding the logRatio for each chromosome
    
    modes = [] ##Modes of log_ratios for each chromosome to take care of the systematic bias
    
    
    ## Making a list of accepted chromosomes before calculating the log ratios
    chrs1 = ['%s' %i for i in range(1, 23)]
    chrs1.append('X')
    chrs1.append('Y')
    chrs2 = ['chr%s' %i for i in range(1, 23)]
    chrs2.append('chrX')
    chrs2.append('chrY')
    accepted_chroms = chrs1 + chrs2
    
    
    for chrom in reference_sample:
        if not chrom in accepted_chroms: continue
        
        ref_sample_values = np.array(reference_sample[chrom].values).astype(float)
        
        mask = ref_sample_values < min_read
        ref_sample_values[mask] = np.nan

        tumor_sample_values = np.array(tumor_sample[chrom].values)

        mask = tumor_sample_values == 0
        tumor_sample_values[mask] = np.nan #np.NINF
        
        log_r = np.log2( tumor_sample_values / ref_sample_values )
        
        chr = Chromosome() #chromosome object containing the start, step, and logRatios
        chr.start = tumor_sample[chrom].start
        chr.step = tumor_sample[chrom].step
        chr.values = log_r
        log_ratios[chrom] = chr 

        ######################
        ## Finding the mode ##
        ######################

        bin_width = 0.05
        minimum = 1000      ##minimum log_ratios
        ## Finding the minimum log_ratios   
        minimum = np.nanmin(log_r)   ## nanmin() also takes care of negative and positive infinity and returns the minimum
        ## Making a histogram and finding the bin with maximum count
        finite_log_ratios = log_r[np.isfinite(log_r)]
        bins = np.floor((finite_log_ratios - minimum) / bin_width ).astype('int')
        del finite_log_ratios

        bin_counts = np.bincount(bins)    ## bincount counts number of occurrences of each value in array of non-negative ints.
        #bin_counts = sig.bincount(bins)    ## bincount counts number of occurrences of each value in array of non-negative ints.
        
        max_bin = np.argmax(bin_counts)
        del bins

        mode = max_bin * bin_width + (bin_width / 2) + minimum  ## The mode of logRatios for a chromosome
        modes.append(mode)
    
    #delete_mitochondrial_chromosome(log_ratios)

    ##Taking care of the systematic bias
    # sys_bias = np.median(modes)
    sys_bias = np.median(modes) if len(modes) != 0 else 0   ##this is for when there is no systematic bias
    sys.stderr.write("Estimated systematic bias: %.3f.\n" %(sys_bias))
    
    for chrom in log_ratios:
        log_ratios[chrom].values = log_ratios[chrom].values - sys_bias
    
    sys.stderr.write("minimum number of reads required at each location to calculate the coverage log ratio: %d.\n" %(min_read))
    return log_ratios     ## Bias free log ratios

def replace_middle_NaNs(filtered_logratios):
    """Replaces the middle NaNs with the average of their boundaries"""
    
    for chrom in filtered_logratios:
        for lst1 in NaN_ranges(filtered_logratios[chrom].values):
            if lst1[0]!=0 and lst1[-1] != len(filtered_logratios[chrom].values)-1 and (lst1[-1] - lst1[0]) < 1000:
                avg = (filtered_logratios[chrom].values[lst1[0]-1] + filtered_logratios[chrom].values[lst1[-1]+1]) / 2
                filtered_logratios[chrom].values[lst1[0]:lst1[-1]+1] = avg
    
    return filtered_logratios

##def delete_mitochondrial_chromosome(log_ratios):
##    if 'chrM' in log_ratios:
##        del log_ratios['chrM']
##        sys.stderr.write("chrM was deleted!\n")
##    if 'chrMT' in log_ratios:
##        del log_ratios['chrMT']
##        sys.stderr.write("chrMT was deleted!\n")
##    if 'M' in log_ratios:
##        del log_ratios['M']
##        sys.stderr.write("Chromosome M was deleted!\n")
##    if 'MT' in log_ratios:
##        del log_ratios['MT']
##        sys.stderr.write("Chromosome MT was deleted!\n")
##    if 'chrNC_007605' in log_ratios:
##        del log_ratios['chrNC_007605']
##        sys.stderr.write("Chromosome NC_007605 was deleted!\n")

def find_consensus_change_points_win_wise(log_ratios, allele_fractions, win_size, logr_thresh, baf_thresh):
    change_points = {}
    
    ##Finding the mean differences for coverage log ratios
    log_ratio_mean_diffs = double_sliding_window_all_chromosomes(log_ratios, win_size)
    
    ##extracting the start and the step from log ratios to be used for B allele fractions
    start = log_ratio_mean_diffs[list(log_ratio_mean_diffs.keys())[0]].start
    step = log_ratio_mean_diffs[list(log_ratio_mean_diffs.keys())[0]].step
    
    ##Finding the mean differences for allele fractions
    mean_filtered_allele_fractions = allele_frac_mean_filt(allele_fractions, win_size, step, start)
    allele_fraction_mean_diffs = find_mean_difference(mean_filtered_allele_fractions, win_size)
    del mean_filtered_allele_fractions

    ##Making wig files out of mean differences of allele fractions and log ratios
    # make_wig_file(log_ratio_mean_diffs, 'log_ratio_mean_diffs_' + str(win_size) + '_' + str(logr_thresh) + '_' + str(baf_thresh) + '.wig')
    # make_wig_file(allele_fraction_mean_diffs, 'allele_fraction_mean_diffs_' + str(win_size) + '_' + str(logr_thresh) + '_' + str(baf_thresh) + '.wig')

    ##Finding the change points for one window size    
    for chrom in log_ratio_mean_diffs:

        logr = np.array(log_ratio_mean_diffs[chrom].values)
        if chrom not in allele_fraction_mean_diffs: continue    ##this was added on 13 JANUARY 2016 because there might be cases where 1 complete chromosome does not have any heterozygous snps!
        
        baf = np.array(allele_fraction_mean_diffs[chrom].values)   

        if len(baf) < len(logr):
            baf = np.concatenate((baf, np.zeros( len(logr) - len(baf))) )
        elif len(baf) > len(logr):
            logr = np.concatenate((logr, np.zeros( len(baf) - len(logr)) ) )

        ## Replacing NaNs with zeros
        logr = np.nan_to_num(logr)
        baf = np.nan_to_num(baf)
        
        compound_scores = (logr / logr_thresh) ** 2 + (baf / baf_thresh) ** 2
        del logr; del baf

        ## Replacing compound_scores below one with zero  (decision boundary is a circle)
        mask = compound_scores < 1
        compound_scores[mask] = 0
        indices = sig.argrelmax(compound_scores)[0]  ## Finds the local maxima (index of local change points with maximum value)   ## To find the single change point

        ##############################################
        ## Fix for when no change point is detected ##
        ##############################################
        if len(indices) == 0:
            sys.stderr.write("Detected no break points for %s with win_size: %s; clogr_threshold: %s, and BAF_threshold: %s.\n" %(chrom, win_size, logr_thresh, baf_thresh))
            chr = Chromosome()
            chr.start = log_ratio_mean_diffs[chrom].start
            chr.step = log_ratio_mean_diffs[chrom].step
            chr.values = []
            change_points[chrom] = chr
            continue
            
        #####################
        ## Fixing the gaps ##
        #####################
        zipped = zip(compound_scores[indices], indices)
        sorted_local_maxima = sorted(zipped, reverse = True)
        mask = np.zeros(len(compound_scores), dtype = np.bool)
        break_points = []
        for max in sorted_local_maxima:
            if mask[max[1]] == False:
                mask[ max[1] - win_size : max[1] + win_size + 1 ] = True
                break_points.append(max[1])

        chr = Chromosome()
        chr.start = log_ratio_mean_diffs[chrom].start
        chr.step = log_ratio_mean_diffs[chrom].step
        chr.values = break_points

        change_points[chrom] = chr
        
    return change_points

def NaN_ranges(lst):
    """Finds the range of consecutive NaNs in a list"""
    idx = list(np.where(np.isnan(lst))[0])
    return np.array_split(idx, list(np.where(np.diff(idx)!= 1)[0]+1))  

def NaN_boundaries(lst):
    """Returns the list of the boundaries of NaNs at the start, centromere, and the end of the chromosome"""
    ranges = []
    for l in NaN_ranges(lst):
        if len(l) == 0: continue
        if l[0] == 0 or l[-1] == len(lst)-1 or (l[-1] - l[0]) >= 1000:
            ranges.append([l[0], l[-1]])
    return ranges

def consensus(change_points_list, window_sizes, log_ratios):
    """finds the consensus change points"""
    ##list containing several change points based on baf and clogr
    consensus_cp = {}
    
    for chrom in change_points_list[0]:
        change_points_indices = [cp[chrom].values for cp in change_points_list]
        ## total: keeps the indices of the change points in the list of values
        ## Populating 'total' with the change points from the smallest win_size
        
        tot_set = set(change_points_indices[0])         

        for i in range(0, len(window_sizes) - 1):
            vicinity = window_sizes[i+1]
            temp = []
            for index in change_points_indices[i+1]: 
                range_index = set(range(index - vicinity, index + vicinity + 1)) ##fixing the exclusiveness of the range() by adding 1
                if len(tot_set.intersection(range_index)) == 0:
                    temp.append(int(index))
            tot_set = tot_set.union(set(temp))    
        del change_points_indices
        total = list(tot_set)
        total.append(0) 
        total.append(len(log_ratios[chrom].values - 1))
        
        for nan_range in NaN_boundaries(log_ratios[chrom].values):
            total.append(nan_range[0])
            total.append(nan_range[1])
        
        chr = Chromosome()
        chr.start = change_points_list[0][chrom].start
        chr.step = change_points_list[0][chrom].step
        chr.values = sorted(set(total))
        consensus_cp[chrom] = chr    
    return consensus_cp    

def find_consensus_only_for_logr_win_wise(log_ratios, win_size, logr_thresh):
    change_points = {}
    log_ratio_mean_diffs = double_sliding_window_all_chromosomes(log_ratios, win_size)
    
    start = log_ratio_mean_diffs[list(log_ratio_mean_diffs.keys())[0]].start
    step = log_ratio_mean_diffs[list(log_ratio_mean_diffs.keys())[0]].step
    
    for chrom in log_ratio_mean_diffs:
        logr = np.array(log_ratio_mean_diffs[chrom].values)
        
        ## Replacing NaNs with zero
        logr = np.nan_to_num(logr)
        
        mask = logr < logr_thresh
        logr[mask] = 0
        indices = sig.argrelmax(logr)[0]
        
        if len(indices) == 0:
            sys.stderr.write("Detected no break points for %s with win_size: %s; clogr_threshold: %s.\n" %(chrom, win_size, logr_thresh))
            chr = Chromosome()
            chr.start = log_ratio_mean_diffs[chrom].start
            chr.step = log_ratio_mean_diffs[chrom].step
            chr.values = []
            change_points[chrom] = chr
            continue
            
        zipped = zip(logr[indices], indices)
        sorted_local_maxima = sorted(zipped, reverse = True)
        mask = np.zeros(len(logr), dtype = np.bool)
        break_points = []
        for max in sorted_local_maxima:
            if mask[max[1]] == False:
                mask[ max[1] - win_size : max[1] + win_size + 1 ] = True
                break_points.append(max[1])

        chr = Chromosome()
        chr.start = log_ratio_mean_diffs[chrom].start
        chr.step = log_ratio_mean_diffs[chrom].step
        chr.values = break_points

        change_points[chrom] = chr
    
    return change_points

def make_seg(log_ratios, allele_fractions, change_points, sample_id, logr_merge_thresh, baf_merge_thresh, to_print, file_name_descriptor):
    """make a seg file containing the segments based on cLOGr and BAF"""
    if to_print == 'True': 
        sys.stderr.write("Printing the results to the standard output...\n")
    else: 
        sys.stderr.write("Printing the results to file named: %s...\n" %(sample_id + file_name_descriptor + '.seg'))
    sys.stderr.write("coverage log ratio merging threshold: %.2f.and B-allele fraction merging threshold: %.2f\n" %(logr_merge_thresh, baf_merge_thresh))
    if to_print == 'False':
        file = open(sample_id + file_name_descriptor + '.seg', 'w')
        file.write("'ID\tchrom\tloc.start\tloc.end\tnum_probes\tseg.mean\tbaf.mean\n")
    else:
        print("'ID\tchrom\tloc.start\tloc.end\tnum_probes\tseg.mean\tbaf.mean")
    
    for chrom in change_points:
        start = int(change_points[chrom].start)
        step  = int(change_points[chrom].step)
        
        cp_indices = change_points[chrom].values
        
        log_ratio_values = log_ratios[chrom].values
        
        allele_fraction_values = allele_fractions[chrom].values #mean_filtered_al[chrom].values
        allele_fraction_positions = allele_fractions[chrom].positions
        
        ii = 0 
        i = 0  
        while i < len(cp_indices) - 2:

            ## Getting rid of the 'empty slice' warning
            mask = np.isfinite(log_ratio_values[cp_indices[i]:cp_indices[i+1]])
            if len(log_ratio_values[cp_indices[i]:cp_indices[i+1]][mask]) == 0:
               clogr_avg1 = np.nan
            else:
               clogr_avg1 =  np.nanmean(log_ratio_values[cp_indices[i]:cp_indices[i+1]])
            
            mask = np.isfinite(log_ratio_values[cp_indices[i+1]:cp_indices[i+2]])
            if len(log_ratio_values[cp_indices[i+1]:cp_indices[i+2]][mask]) == 0:
                clogr_avg2 = np.nan
            else:
                clogr_avg2 = np.nanmean(log_ratio_values[cp_indices[i+1]:cp_indices[i+2]])
            
            
            first_first = ii    
            while allele_fraction_positions[ii] < cp_indices[i] * step and ii < len(allele_fraction_positions) - 1:  ##start
                ii += 1
            j = ii 
            while allele_fraction_positions[j] <= cp_indices[i+1] * step and j < len(allele_fraction_positions) - 1:  ##end
                j += 1
            
            ## Getting rid of the 'empty slice' warning
            mask = np.isfinite(allele_fraction_values[ii:j])
            if len(allele_fraction_values[ii:j][mask]) == 0:
                baf_avg1 = np.nan
            else:
                baf_avg1 = np.nanmean(allele_fraction_values[ii:j])

                
            first_end = j
            ii = j   ##end of previous is the start of the current
            while allele_fraction_positions[ii] < cp_indices[i+1] * step and ii < len(allele_fraction_positions) - 1:
                ii += 1
            j = ii 
            while allele_fraction_positions[j] <= cp_indices[i+2] * step and j < len(allele_fraction_positions) - 1:
                j += 1
            
            ## Getting rid of the 'empty slice' warning
            mask = np.isfinite(allele_fraction_values[ii:j])
            if len(allele_fraction_values[ii:j][mask]) == 0:
                baf_avg2 = np.nan
            else:
                baf_avg2 = np.nanmean(allele_fraction_values[ii:j])
            
            ii = first_end 
            
            ##merging
            if str(clogr_avg1) == 'nan':
                if to_print == 'False':
                    file.write("%s\t%s\t%d\t%d\t%d\t%.4f\t%.4f\n" %(sample_id, chrom, (start +(cp_indices[i] * step)), (start + (cp_indices[i+1] * step)), int(((start + (cp_indices[i+1] * step)) - (start +(cp_indices[i] * step)))/step), clogr_avg1, baf_avg1))
                else:
                    print("%s\t%s\t%d\t%d\t%d\t%.4f\t%.4f" %(sample_id, chrom, (start +(cp_indices[i] * step)), (start + (cp_indices[i+1] * step)), int(((start + (cp_indices[i+1] * step)) - (start +(cp_indices[i] * step)))/step), clogr_avg1, baf_avg1))
                i += 1
            elif str(clogr_avg2) == 'nan':
                if to_print == 'False':
                    file.write("%s\t%s\t%d\t%d\t%d\t%.4f\t%.4f\n" %(sample_id, chrom, (start +(cp_indices[i] * step)), (start + (cp_indices[i+1] * step)), int(((start + (cp_indices[i+1] * step)) - (start +(cp_indices[i] * step)))/step), clogr_avg1, baf_avg1))
                else:
                    print("%s\t%s\t%d\t%d\t%d\t%.4f\t%.4f" %(sample_id, chrom, (start +(cp_indices[i] * step)), (start + (cp_indices[i+1] * step)), int(((start + (cp_indices[i+1] * step)) - (start +(cp_indices[i] * step)))/step), clogr_avg1, baf_avg1))
                i += 1
            elif abs(clogr_avg2 - clogr_avg1) > logr_merge_thresh or abs(baf_avg2 - baf_avg1) > baf_merge_thresh:
                if to_print == 'False':
                    file.write("%s\t%s\t%d\t%d\t%d\t%.4f\t%.4f\n" %(sample_id, chrom, (start +(cp_indices[i] * step)), (start + (cp_indices[i+1] * step)), int(((start + (cp_indices[i+1] * step)) - (start +(cp_indices[i] * step)))/step), clogr_avg1, baf_avg1))
                else:
                    print("%s\t%s\t%d\t%d\t%d\t%.4f\t%.4f" %(sample_id, chrom, (start +(cp_indices[i] * step)), (start + (cp_indices[i+1] * step)), int(((start + (cp_indices[i+1] * step)) - (start +(cp_indices[i] * step)))/step), clogr_avg1, baf_avg1))
                i += 1
            else:
                cp_indices.remove(cp_indices[i+1])   ##If non satisfies the conditions, they are merged
                ii = first_first
        
        if to_print == 'False':
            file.write("%s\t%s\t%d\t%d\t%d\t%.4f\t%.4f\n" %(sample_id, chrom,(start + (cp_indices[i] * step)), (start + (cp_indices[i+1] * step)), int(((start + (cp_indices[i+1] * step)) - (start +(cp_indices[i] * step)))/step), clogr_avg2, baf_avg2))
        else:
            print("%s\t%s\t%d\t%d\t%d\t%.4f\t%.4f" %(sample_id, chrom,(start + (cp_indices[i] * step)), (start + (cp_indices[i+1] * step)), int(((start + (cp_indices[i+1] * step)) - (start +(cp_indices[i] * step)))/step), clogr_avg2, baf_avg2))
        ##end of merging 
    if to_print == 'False': file.close()

def print_seg_only_for_coverage(log_ratios, change_points, sample_id, logr_merge_thresh, to_print, file_name_descriptor):
    if to_print == 'True': 
        sys.stderr.write("Printing the results to the standard output...\n")
    else: 
        sys.stderr.write("Printing the results to file named: %s...\n" %(sample_id + file_name_descriptor + '.seg'))
    sys.stderr.write("coverage log ratio merging threshold: %.2f.\n" %(logr_merge_thresh))
    if to_print == 'False':
        file = open(sample_id + file_name_descriptor + '.seg', 'w')
        file.write("'ID\tchrom\tloc.start\tloc.end\tnum_probes\tseg.mean\n")
    else:
        print("'ID\tchrom\tloc.start\tloc.end\tnum_probes\tseg.mean")
    
    for chrom in change_points:
        start = int(change_points[chrom].start)
        step  = int(change_points[chrom].step)
        
        cp_indices = change_points[chrom].values
        
        log_ratio_values = log_ratios[chrom].values

        i = 0  
        while i < len(cp_indices) - 2:
            ## Getting rid of the 'empty slice' warning
            mask = np.isfinite(log_ratio_values[cp_indices[i]:cp_indices[i+1]])
            if len(log_ratio_values[cp_indices[i]:cp_indices[i+1]][mask]) == 0:
               clogr_avg1 = np.nan
            else:
               clogr_avg1 =  np.nanmean(log_ratio_values[cp_indices[i]:cp_indices[i+1]])
            
            mask = np.isfinite(log_ratio_values[cp_indices[i+1]:cp_indices[i+2]])
            if len(log_ratio_values[cp_indices[i+1]:cp_indices[i+2]][mask]) == 0:
                clogr_avg2 = np.nan
            else:
                clogr_avg2 = np.nanmean(log_ratio_values[cp_indices[i+1]:cp_indices[i+2]])

            if str(clogr_avg1) == 'nan':
                if to_print == 'False':
                    file.write("%s\t%s\t%d\t%d\t%d\t%.4f\n" %(sample_id, chrom, (start +(cp_indices[i] * step)), (start + (cp_indices[i+1] * step)), int(((start + (cp_indices[i+1] * step)) - (start +(cp_indices[i] * step)))/step), clogr_avg1))
                else:
                    print("%s\t%s\t%d\t%d\t%d\t%.4f" %(sample_id, chrom, (start +(cp_indices[i] * step)), (start + (cp_indices[i+1] * step)), int(((start + (cp_indices[i+1] * step)) - (start +(cp_indices[i] * step)))/step), clogr_avg1))
                i += 1
            elif str(clogr_avg2) == 'nan':
                if to_print == 'False':
                    file.write("%s\t%s\t%d\t%d\t%d\t%.4f\n" %(sample_id, chrom, (start +(cp_indices[i] * step)), (start + (cp_indices[i+1] * step)), int(((start + (cp_indices[i+1] * step)) - (start +(cp_indices[i] * step)))/step), clogr_avg1))
                else:
                    print("%s\t%s\t%d\t%d\t%d\t%.4f" %(sample_id, chrom, (start +(cp_indices[i] * step)), (start + (cp_indices[i+1] * step)), int(((start + (cp_indices[i+1] * step)) - (start +(cp_indices[i] * step)))/step), clogr_avg1))
                i += 1
            elif abs(clogr_avg2 - clogr_avg1) > logr_merge_thresh:
                if to_print == 'False':
                    file.write("%s\t%s\t%d\t%d\t%d\t%.4f\n" %(sample_id, chrom, (start +(cp_indices[i] * step)), (start + (cp_indices[i+1] * step)), int(((start + (cp_indices[i+1] * step)) - (start +(cp_indices[i] * step)))/step), clogr_avg1))
                else:
                    print("%s\t%s\t%d\t%d\t%d\t%.4f" %(sample_id, chrom, (start +(cp_indices[i] * step)), (start + (cp_indices[i+1] * step)), int(((start + (cp_indices[i+1] * step)) - (start +(cp_indices[i] * step)))/step), clogr_avg1))
                i += 1
            else:
                cp_indices.remove(cp_indices[i+1])   ##If non satisfies the conditions, they are merged
        if to_print == 'False':
            file.write("%s\t%s\t%d\t%d\t%d\t%.4f\n" %(sample_id, chrom, (start + (cp_indices[i] * step)), (start + (cp_indices[i+1] * step)), int(((start + (cp_indices[i+1] * step)) - (start +(cp_indices[i] * step)))/step), clogr_avg2))
        else:
            print("%s\t%s\t%d\t%d\t%d\t%.4f" %(sample_id, chrom, (start + (cp_indices[i] * step)), (start + (cp_indices[i+1] * step)), int(((start + (cp_indices[i+1] * step)) - (start +(cp_indices[i] * step)))/step), clogr_avg2))
    if to_print == 'False': file.close()

def find_allele_fraction_segment_mean(allele_fraction_positions, allele_fraction_values, start, end):
    """"""
    i = 0 
    while allele_fraction_positions[i] < start:
        i += 1
    j = i 
    while allele_fraction_positions[j] <= end:
        j += 1
    return np.nanmean(allele_fraction_values[i:j])

def make_seg_only_logr(seg_file):
    file = open(seg_file, 'r')
    file2 = open('clogr_' + seg_file, 'w')
    header = file.readline().split('\t')
    s = header[0] + '\t' + header[1] + '\t' + header[2] + '\t' + header[3] + '\t' + header[4]
    file2.write(s + '\n')
    for line in file:
        s = line.split('\t')
        s2 = s[0] + '\t' + s[1] + '\t' + s[2] + '\t' + s[3] + '\t' + s[4]
        file2.write(s2 + '\n')
    file.close()
    file2.close()   

def make_seg_only_baf(seg_file):
    file = open(seg_file, 'r')
    file2 = open('baf_' + seg_file, 'w')
    header = file.readline().split('\t')
    s = header[0] + '\t' + header[1] + '\t' + header[2] + '\t' + header[3] + '\t' + header[5].strip()
    file2.write(s + '\n')
    for line in file:
        s = line.split('\t')
        s2 = s[0] + '\t' + s[1] + '\t' + s[2] + '\t' + s[3] + '\t' + s[5].strip()
        file2.write(s2 + '\n')
    file.close()
    file2.close()

def segmentum_with_baf(min_read, sample_file, ref_file, snp_file, win_size, clogr_threshold, baf_thresh, logr_merge_thresh, baf_merge_thresh, z_score, to_print, BAF_to_WIG):
    st = time.time()
    log_ratios = calculate_logratios(read_wig(sample_file), read_wig(ref_file), min_read)

    ## Deleting ChrY for the moment...
    ## This is because for the moment we don't get any signal from chrY
    ## because of no heterozygous snp is detected.
    log_ratios.pop('chrY', None)
    filtered_log_ratios = filter_outliers(log_ratios, 5)
    filtered_log_ratios = replace_middle_NaNs(filtered_log_ratios)
    del log_ratios
    sample_file_name = os.path.basename(os.path.splitext(sample_file)[0])
    if 'wig' in sample_file_name: sample_file_name = sample_file_name.replace('.wig', '')
    allele_fractions = extract_sample(snp_file, sample_file_name)
    smoothed_allele_fractions = smooth_allele_fractions(allele_fractions, 9)
    del allele_fractions
    
    ## extracting the allele fractions to check the found regions with cnLOH
    if BAF_to_WIG:
        baf_file_name= sample_file_name + '_BAFs.wig'
        make_wig_file2(smoothed_allele_fractions, baf_file_name)
    
    win_sizes, logr_thresholds, baf_thresholds = estimate_threshold(win_size, clogr_threshold, 
            baf_thresh, logr_merge_thresh, baf_merge_thresh, z_score)
    
    change_points_list = [find_consensus_change_points_win_wise(filtered_log_ratios, smoothed_allele_fractions,
           win_sizes[i], logr_thresholds[i], baf_thresholds[i]) for i in range(0,len(win_sizes))]
    
    consensus_change_points = consensus(change_points_list, win_sizes, filtered_log_ratios)
    file_name_descriptor = "_%s_%s_%s" %(win_size, clogr_threshold, baf_thresh)
    
    make_seg(filtered_log_ratios, smoothed_allele_fractions, consensus_change_points, 
        sample_file_name, logr_merge_thresh, baf_merge_thresh, to_print, file_name_descriptor)
    ## make_seg_only_logr(sample_file_name + '.seg')
    ## make_seg_only_baf(sample_file_name + '.seg')
    m, s = divmod(time.time() - st, 60)
    h, m = divmod(m, 60)
    sys.stderr.write("Total time: %d hours, %02d minutes and %02d seconds.\n" %(h, m, s))

def segmentum_without_baf(min_read, sample_file, ref_file, win_size, clogr_threshold, logr_merge_thresh, z_score, to_print):
    st = time.time()
    log_ratios = calculate_logratios(read_wig(sample_file), read_wig(ref_file), min_read)
    filtered_log_ratios = filter_outliers(log_ratios, 5)
    filtered_log_ratios = replace_middle_NaNs(filtered_log_ratios)
    del log_ratios
    sample_file_name = os.path.basename(os.path.splitext(sample_file)[0])
    if 'wig' in sample_file_name: sample_file_name = sample_file_name.replace('.wig', '')
    file_name_descriptor = "_%s_%s" %(win_size, clogr_threshold)
    win_sizes, logr_thresholds = estimate_threshold_only_for_coverage(win_size, clogr_threshold, logr_merge_thresh, z_score)
    change_points_list = [find_consensus_only_for_logr_win_wise(filtered_log_ratios, win_sizes[i], logr_thresholds[i]) for i in range(0, len(win_sizes))]
    consensus_change_points = consensus(change_points_list, win_sizes, filtered_log_ratios)
    print_seg_only_for_coverage(filtered_log_ratios, consensus_change_points, sample_file_name, logr_merge_thresh, to_print, file_name_descriptor)
    m, s = divmod(time.time() - st, 60)
    h, m = divmod(m, 60)
    sys.stderr.write("Total time: %d hours, %02d minutes and %02d seconds.\n" %(h, m, s))


if __name__ == '__main__':
    args = docopt.docopt(__doc__)
    if args['extract'] and args['read'] and args['depth']:
        wsize = int(args['<window_size>'])
        step = wsize / 2
        mode = 'A'
        coverage_tiled(args['<BAM_file>'], wsize, quality=int(args['--quality']), mode=mode, step=step)
    elif args['calculate'] and args['BAF']:
        options = Options(args)
        bam_files = [args['<tumor>'], args['<normal>']]
        calculate_BAF(bam_files, args['<genome_fasta>'], args['<SNP_position>'], options)
    elif args['analyze'] and args['with'] and args['BAF']:
        segmentum_with_baf(int(args['--min_read']), args['<tumor>'], args['<normal>'], args['<BAF_file>'], int(args['<window_size>']), float(args['<clogr_threshold>']), float(args['<BAF_threshold>']), float(args['--logr_merge']), float(args['--baf_merge']), 5, args['--print'], args['--BAF'])
    elif args['analyze'] and args['without'] and args['BAF']:
        segmentum_without_baf(int(args['--min_read']), args['<tumor>'], args['<normal>'], int(args['<window_size>']),  float(args['<clogr_threshold>']), float(args['--logr_merge']), 5, args['--print'])
    elif args['find'] and args['recurrent'] and ['cnLOHs']:
        chrom_lengths = find_chromosome_length(args['<seg_files>'])
        clogr_thresh = float(args['--clogr_thresh'])
        baf_thresh = float(args['--baf_thresh'])
        loh_regions_chromosome_wise = find_cnloh(args['<seg_files>'], clogr_thresh, baf_thresh)
        find_recurrent(loh_regions_chromosome_wise, chrom_lengths)
    elif args['simulate']:
        normal_contamination = 1 - float(args['--tumor_purity'])
        output_prefix = args['--output_prefix']
        read_length = int(args['--read_length'])
        file_path = args['<normal>']
        simulate(normal_contamination, output_prefix, read_length, file_path)
