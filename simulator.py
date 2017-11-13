#!/usr/bin/env python3
"""
A tool for simulating normal and tumour read depth and B-allele fraction of Heterozygous SNPs.

Usage:
  simulator <normal_file> [-c N] [-p N] [-r N]

Options:
  -h --help         Show this screen.
  -c --normal_contamination=N    percentage of contamination by normal tissue [default: 0.3]
  -p --output_prefix=N    prefix to be assigned to the simulated files [default: simulated]
  -r --read_length=N    read length to be considered for simulation [default: 150]

Author: Ebrahim Afyounian <ebrahim.afyounian@uta.fi>
"""
import time, docopt, sys, re, gzip, subprocess
import numpy as np
import scipy.signal as sig
from scipy import stats
#import matplotlib
#matplotlib.use('Agg')
#from mpl_toolkits.mplot3d import Axes3D
#from matplotlib import cm
#from matplotlib.ticker import LinearLocator, FormatStrFormatter
#import matplotlib.pyplot as plt

class Chromosome:
    """class to keep track of start, step, and (coverage OR logRatios) data from a wig file"""
    def __init__(self):
        self.start = None 
        self.step = None
        self.values = [] 

class AlleleFraction:
    """"""
    def __init__(self):
        self.values = [] 
        self.positions = [] 

def open_file(in_path):
    if in_path.endswith('.gz'):
        file = gzip.open(in_path, 'r')
    else:
        file = open(in_path, 'r')
    return file 

def read_wig(sampleFileName):
    """reads a wig file and returns a dictionary out of it"""
    if sampleFileName.endswith('.gz'):
        file = subprocess.Popen('gunzip -c %s' % sampleFileName, stdout=subprocess.PIPE, shell=True).stdout
    else:
        file = open(sampleFileName, 'r')

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
            if float(line) == 0:
                values.append(np.nan)
            else:    
                values.append(float(line))
    file.close()   
    return chromosomes

def make_wig(sample, out_path):
    file = open(out_path, 'w')
    file.write('track graphType=bar viewLimits=0:1100 windowingFunction=none\n')    
    for chrom in sample:
        start = sample[chrom].start
        step = sample[chrom].step
        values = sample[chrom].values
        
        file.write("fixedStep chrom=%s start=%s step=%s\n" %(chrom, start, step))

        for value in values:
            if str(value) == 'nan':
                file.write('NaN' + '\n')
            else:
                file.write(str(value) + '\n')
    file.close()

def flip_coin():
    return np.random.choice(['head', 'tail'], 1)[0]

def calculate_BAF(a,b):
    c = 1 if (a+b) == 0 else a+b
    if np.random.choice([0,1]) == 1:
        return float(a) / c
    else:
        return float(b) / c

def find_max_read_depth(file):
    max = 0; sum = 0; count = 0
    for line in file:
        line = line.decode('utf8')
        if line.startswith("fixedStep") or line.startswith("variableStep"): continue
        if int(line) != 0:
            sum += int(line)
            count += 1
    mean = float(sum) / count
    max = 2 * mean
    file.close()
    return int(max)

def populate_affine_mat(affine_matrix, file, telo_centromeres):

    dim = affine_matrix.shape[0]
    for line in file:
        line = line.decode('utf8')
        if line[0] in 'fv': 
            chrom = re.search('chrom=(\S+)', line).group(1)
            step = int(re.search('step=(\d+)', line).group(1))
            line_count = 0
            next_value = -1
            continue
         
        pre_value = next_value
        next_value = int(line.strip())
        line_count += 1
        if chrom in ('chrX', 'chrY', 'chrM', 'chrMT'): continue
        if pre_value >= dim or next_value >= dim: continue
        in_range = False

        for rangee in telo_centromeres[chrom]:
            if rangee[0] - 500000 <= line_count * step <= rangee[1] + 500000:
                in_range = True
                break
        
        if not in_range and pre_value >= 0: 
            affine_matrix[pre_value, next_value] += 1
    file.close()
    return affine_matrix

def draw_surface(pop_affine_matrix, max_num_read):
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    X = np.arange(0, max_num_read - 1, 1)
    Y = np.arange(0, max_num_read - 1, 1)
    Z = pop_affine_matrix
    # Z[0,0] = 0   ##should be removed after I take care of the centromere and telomere zeros
    X, Y = np.meshgrid(X, Y)
    surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,
        linewidth=0, antialiased=False)
    fig.colorbar(surf)
    # plt.savefig('affin_landscape.png')
    plt.savefig('affin_landscape.svg')

def delete_mitochondrial_chromosome(logRatioDict):
    if 'chrM' in logRatioDict.keys():
        del logRatioDict['chrM']
    if 'chrMT' in logRatioDict.keys():
        del logRatioDict['chrMT']

def telo_centro(coverage):
    telo_centromeres = {}
    for chrom in coverage:
            step = int(coverage[chrom].step)
            values = coverage[chrom].values
            indices = list(np.where(np.isnan(values))[0])
            nan_ranges = np.array_split(indices, list(np.where(np.diff(indices)!= 1)[0]+1))
            nan_ranges_greater_than_50kb = []
            for l in nan_ranges:
                if l[0] == 0 or l[-1] == len(values)-1 or (l[-1] - l[0]) >= 50:
                    nan_ranges_greater_than_50kb.append((l[0] * step, l[-1] * step))   ## tuple of (start, end)
            telo_centromeres[chrom] = nan_ranges_greater_than_50kb    
    return telo_centromeres

def extract_chromosome_sizes(coverage):
    chrom_sizes = {}
    for chrom in coverage:
        step = int(coverage[chrom].step)
        chrom_sizes[chrom] = (len(coverage[chrom].values) -1) * step
    return chrom_sizes

def smirnov_transform(affine_matrix, max_num_read):
    """performs an inverse transform sampling"""
    rv_discretes = []
    xk = np.arange(int(max_num_read))
    for i in range(0, (int(max_num_read))):
        pk = affine_matrix[:,i]
        rv_discretes.append(stats.rv_discrete(name='custm', values=(xk, pk)))    
    sys.stderr.write('Smirnov transform done.\n')
    return rv_discretes

def simulate_RD(chrom_sizes, step, rv_discretes):
    pre_normal = {}
    for chrom in chrom_sizes:
        chrom_size = int(chrom_sizes[chrom]) 
        start_value = 0
        values = [start_value]
        pre_value = start_value
        count = 0 
        while count <= chrom_size / step:
            next_value = rv_discretes[int(pre_value)].rvs(size=1)
            values.append(next_value)
            pre_value = next_value
            count += 1
        chr = Chromosome()
        chr.start = 0
        chr.step = step
        chr.values = values
        
        pre_normal[chrom] = chr
    return pre_normal            

def simulate(normal_contamination, output_prefix, read_length, file_path):

    st = time.time()
    file = open_file(file_path)
    max_num_read = find_max_read_depth(file)
    
    ## Populating the affinity matrix
    affine_matrix = np.zeros((max_num_read, max_num_read), dtype='int32')
    # coverage = read_wig(args['<normal_file>'])
    coverage = read_wig(file_path)
    delete_mitochondrial_chromosome(coverage)
    telo_centromeres = telo_centro(coverage)
    # sys.stderr.write('telomeres and centromeres were found.\n')

    chrom_sizes = extract_chromosome_sizes(coverage)
    # sys.stderr.write('chromosome sizes were found.\n')
    del coverage
    
    # file = open_file(args['<normal_file>'])
    file = open_file(file_path)
    pop_affine_matrix =  populate_affine_mat(affine_matrix, file, telo_centromeres)
    # sys.stderr.write('affine matrix: populated.\n')
    
    ## Visualization
    # draw_surface(pop_affine_matrix, max_num_read)
    ## Normalizing the affinity matrix
    norm_affine_matrix = pop_affine_matrix / np.sum(pop_affine_matrix, axis=0).astype('float')
    # sys.stderr.write('affine matrix: normalized\n')

    ## Using Smirnov transform (Inverse transform sampling) 
    rv_discretes = smirnov_transform(norm_affine_matrix, max_num_read)
    
    ## Extracting the resolution of sequencing
    # file = open_file(args['<normal_file>'])
    file = open_file(file_path)
    step = int(re.search('step=(\d+)', file.readline().decode('utf8')).group(1))
    file.close()
    
    ## Simulation of the normal genome
    sys.stderr.write('Started simulating the pre_normal genome...\n')



    pre_normal = simulate_RD(chrom_sizes, step, rv_discretes)
    sys.stderr.write('Done with simulating the pre_normal genome.\n')
    
    ## Creating the simulated affinity (without noise)
    win_size = 3
    for chrom in pre_normal:
        pre_normal[chrom].values = sig.medfilt(pre_normal[chrom].values, win_size)
    
    #########################    
    ## Creating the normal ##
    #########################

    normal = {}     
    for chrom in pre_normal:        
        chr = Chromosome()
        chr.step = pre_normal[chrom].step
        chr.start = pre_normal[chrom].start
        chr.values = np.random.poisson(pre_normal[chrom].values.astype(int))
        normal[chrom] = chr

    #########################
    ## Creating the tumour ##   
    #########################

    tumour = {}
    allele_fractions = {}
    ground_truth = {}

    ##normal_contamination = 0.3     #assuming 30% of normal cell contamination
    sys.stderr.write('%.2f%% of normal contamination will be used for simulation.\n' %(normal_contamination*100) )
    #read_length = 150
    mean_snp_distance = 1500   # assuming the there are one snp per 1500 bp in the human genome


    for chrom in pre_normal:
        chrom_size = len(pre_normal[chrom].values)

        start = np.random.choice(chrom_size - 10000 , 40, replace = False)
        start.sort()

        lens = []

        lens = np.concatenate((np.random.choice(range(5, 100), 10), lens))
        lens = np.concatenate((np.random.choice(range(101, 500), 10), lens))
        lens = np.concatenate((np.random.choice(range(501, 1000), 10), lens))
        lens = np.concatenate((np.random.choice(range(1001,10000), 10), lens))
        lens = lens.astype(int)
        
        ## added this on 26.09.2014
        lens = np.random.permutation(lens)   ## permutes the lens so that they are randomly distributed throughout the chromosome
        ###########################
        
        end = start + lens

        cn_alterations = np.random.choice([1, -1], 40)

        cn_track_a = np.repeat(1.0, chrom_size)
        cn_track_b = np.repeat(1.0, chrom_size)

        current_chrom_depth = normal[chrom].values

        ###########################################################
        ## Fix for the start of the chr which read depth is zero ##
        ###########################################################
        i = 0
        while current_chrom_depth[i] < 50:    ## < 50 because the minimum read is 50
            i += 1
        start = np.concatenate(([0], start))
        end = np.concatenate(([i], end))
        cn_alterations = np.concatenate(([1], cn_alterations))  ## 1 bcz of line 166 and 167 np.diff
        ###########################################################

        for i in range(0,len(start)):
            if flip_coin() == 'head':   #if 'head' change track a
                if cn_alterations[i] == -1:        
                    cn_track_a[start[i]:end[i]+1][cn_track_a[start[i]:end[i]+1] != 0] += cn_alterations[i]
                else:
                    cn_track_a[start[i]:end[i]+1] += cn_alterations[i]
                
            else:  # if 'tail' change track b    
                if cn_alterations[i] == -1:
                    cn_track_b[start[i]:end[i]+1][cn_track_b[start[i]:end[i]+1] != 0] += cn_alterations[i]
                else:
                    cn_track_b[start[i]:end[i]+1] += cn_alterations[i]

        cnt = cn_track_a + cn_track_b            

        ###############################
        ## Creating the ground truth ##
        ###############################
        
        breaks_in_a = np.sort(np.where(np.diff(cn_track_a))[0])
        breaks_in_b = np.sort(np.where(np.diff(cn_track_b))[0])

        all_breaks = sorted(set(np.concatenate((breaks_in_a, breaks_in_b))))
        
        b1 = np.concatenate(([0], all_breaks))
        b2 = np.concatenate((all_breaks, [chrom_size - 1]))
        segs = zip(b1, b2)
        
        final_copy_numbers = cnt[b2]
        
        segs_and_copy_numbers = zip(segs, final_copy_numbers)
        
        del b1, b2
        ground_truth[chrom] = segs_and_copy_numbers
        del segs, segs_and_copy_numbers
        
        ####################################
        ## creating the tumour read depth ## 
        ####################################
        
        ## introducing the normal cell contamination to the copy number track
        cnt = cnt * (1 - normal_contamination) + (2 * normal_contamination)

        chr = Chromosome()
        chr.start = pre_normal[chrom].start
        chr.step = pre_normal[chrom].step
        chr.values = (np.array(pre_normal[chrom].values) * cnt) / 2
        
        ###################################################
        ## Adding the Poisson noise to the tumour values ##
        ###################################################
        
        for index, value in enumerate(chr.values):
            if str(value) == str(np.nan):
                chr.values[index] = np.nan
            else:    
                chr.values[index] = int(np.random.poisson(int(value), 1))

        tumour[chrom] = chr
        del cnt

        #############################
        ## creating the tumour BAF ## 
        #############################

        step = int(normal[chrom].step)

        ## distributing the snps throughout the chromosome   ##  '- step' is to prevent index out of bound error when rounding
        positions = np.random.choice( (chrom_size - 100) * step, ((chrom_size * step) / mean_snp_distance) - mean_snp_distance, replace=False)  ##chrom_size -100 to fix index n is out of bounds for axis 1 with size n
        positions = sorted(positions)

        ## Calculating the expected number of reads that overlap the snp

        ## Extracting the parameter n in the binomial
        
        n = []    ## Parameter n in the binomial distribution
        p = []    ## Parameter p in the binomial distribution

        ## the index of the nearest datapoint in chr_depth
        ### Vectorization the above code
        indices = np.around(np.array(positions) / step).astype(int)
        indices = indices[:len(indices)-1] ## changed indices to indices[:-1] bcz of:  index n is out of bounds for axis 1 with size n
        number_of_reads = current_chrom_depth[indices].astype(float)  
        n =  read_length * ( number_of_reads / (read_length + step) )
        A = (normal_contamination * 1) + ((1 - normal_contamination) * cn_track_a[indices])
        B = (normal_contamination * 1) + ((1 - normal_contamination) * cn_track_b[indices])
        vfunc = np.vectorize(calculate_BAF)
        p = vfunc(A,B)
        ### End of vectorization
        
        n = np.array(n).astype(int)    
        p = np.array(p)
        baf = np.full(n.size, np.nan)
        mask = n >= 30   # setting to nan those snps which their coverage is less that 30
        baf[mask] = np.random.binomial(n[mask], p[mask]).astype(float) / n[mask]
        
        allele_fraction = AlleleFraction()
        allele_fraction.positions = positions
        allele_fraction.values = baf
        
        allele_fractions[chrom] = allele_fraction

    #############################
    ## Creating the wig files  ##
    #############################

    make_wig(normal, output_prefix + '_normal_rd_' +  str(np.round(normal_contamination, 2)) + '_contamination.wig')
    make_wig(tumour, output_prefix + '_tumor_rd_' +  str(np.round(normal_contamination, 2)) + '_contamination.wig')

    ########################################
    ## Saving the B allele fraction track ##
    ########################################

    file = open(output_prefix + '_BAF' + str(np.round(normal_contamination, 2)) + '_contamination.wig', 'w')
    file.write('track graphType=points viewLimits=0:1 windowingFunction=none\n')
    for chrom in allele_fractions:    
        positions = allele_fractions[chrom].positions
        baf = allele_fractions[chrom].values
        
        file.write('variableStep chrom=' + str(chrom) + '\n')
        positions = positions[: len(positions)-1]  ## changed indices to indices[:-1] bcz of:  index n is out of bounds for axis 1 with size n
        for index, value in enumerate(positions):
            baf_temp = baf[index] if str(baf[index]) != 'nan' else 'NaN'
            file.write(str(positions[index]) + ' ' + str(baf_temp) + '\n')
    file.close()

    ############################
    ## Creating the .tsv file ##
    ############################

    file = open(output_prefix + '_BAF' + str(np.round(normal_contamination, 2)) + '_contamination.tsv', 'w')
    
    col_name = output_prefix + '_tumor_rd_' +  str(np.round(normal_contamination, 2)) + '_contamination.wig'
    file.write('CHROM\tPOSITION\t%s\n' %(col_name))
    for chrom in allele_fractions:
        positions = allele_fractions[chrom].positions
        baf = allele_fractions[chrom].values
        
        for index, value in enumerate(baf):
            if str(value) != 'nan':
                file.write('%s\t%s\t%s\n' %(chrom, positions[index], round(value, 3) ) )
    file.close()

    #############################
    ## Saving the ground truth ##
    #############################
    gt_name  = output_prefix + 'ground_truth_' + str(np.round(normal_contamination, 2)) + '_contamination.seg'
    file = open(gt_name, 'w')
    file.write('\'ID\tchrom\tloc.start\tloc.step\tseg.mean\n')
    for chrom in ground_truth:
        for seg in ground_truth[chrom]:
            file.write('%s\t%s\t%s\t%s\t%s\n' %( 'ground_truth', chrom, int(seg[0][0]) * int(step), int(seg[0][1]) * int(step),  str(seg[1]) ))
    file.close()

    sys.stderr.write('elapsed time: %s\n' %(time.time() - st))

if __name__ == '__main__':
    args = docopt.docopt(__doc__)
    normal_contamination = float(args['--normal_contamination'])
    output_prefix = args['--output_prefix']
    read_length = int(args['--read_length'])
    file_path = args['<normal_file>']
    simulate(normal_contamination, output_prefix, read_length, file_path)
