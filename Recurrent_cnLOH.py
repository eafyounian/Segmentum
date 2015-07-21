#!/bin/env python
"""
A tool for extracting recurrent copy neutral LOH regions across one or many samples.

Usage:
  Recurrent_cnLOH.py <seg_files>... [-c N] [-b N] 

Options:
  -h --help         Show this screen.
  -c --clogr_thresh=N   Coverage logratio must be below this threshold to call a copy neutral LOH region [default: 0.1]
  -b --baf_thresh=N     B-allele fraction must be below this threshold to call a copy neutral LOH region [default: 0.15]   
  
Author: Ebrahim Afyounian <ebrahim.afyounian@uta.fi>
"""

import docopt, itertools, sys
from operator import itemgetter
from quicksect import IntervalNode

def find_chromosome_length(seg_files):
    """finds the chromosome length from the .seg file and returns a dictionary of lengths"""
    chrom_lengths = {}
    chrom = {}
    
    seg_file = open(seg_files[0], 'r')
    header = seg_file.readline()

    for line in seg_file:
        key = line.strip().split('\t')[1]
        seg_end = int(line.strip().split('\t')[3])
        if key not in chrom:
            chrom[key] = [seg_end]
        elif key in chrom:
            chrom[key].append(seg_end)
    
    for key in chrom:
        chrom[key] = sorted(chrom[key])
        chrom_lengths[key] = chrom[key][-1] 
    
    return chrom_lengths

def find_overlaps(start, end, tree):
    """returns overlapping intervals from a interval tree"""
    out = []
    tree.intersect( start, end, lambda x: out.append(x) )
    return [ (x.start, x.end) for x in out ]     
    
def find_cnloh(seg_files, clogr_thresh, baf_thresh):
    """finds regions with loss of heterozygozity in from several .seg files and returns the aggregated loh regions"""
    loh_regions_chromosome_wise = {}
    loh_regions = []
    
    for seg_file in seg_files:
        file = open(seg_file, 'r')
        header = file.readline().split('\t')
        for line in file:
            items = line.strip().split('\t')
            if abs(float(items[4])) <= clogr_thresh and float(items[5]) >= (0.5 - baf_thresh):
                loh_regions.append(items[0:4])
        file.close()
    
    loh_regions = sorted(loh_regions, key=itemgetter(1))
    
    for region in loh_regions:
        if region[1] not in loh_regions_chromosome_wise:
            loh_regions_chromosome_wise[region[1]] = [(int(region[2]),int(region[3]))]
        elif region[1] in loh_regions_chromosome_wise:
            loh_regions_chromosome_wise[region[1]].append((int(region[2]),int(region[3])))

    for chrom in loh_regions_chromosome_wise:
        loh_regions_chromosome_wise[chrom] = sorted(loh_regions_chromosome_wise[chrom], key=itemgetter(0))
    
    return loh_regions_chromosome_wise    

def main(loh_regions_chromosome_wise, chrom_lengths):

    print 'chrom\tseg.start\tseg.end\toccurrence_number' 
    for chrom in loh_regions_chromosome_wise:
        regions = loh_regions_chromosome_wise[chrom]

        ## creating the interval tree
        start, end = regions[0]
        tree = IntervalNode( start, end )
        for start, end in regions[1:]:
            tree = tree.insert( start, end )

        ## generating the query 
        seg_boundary =  set(itertools.chain(*regions))
        seg_boundary.add(0)
        seg_boundary.add(chrom_lengths[chrom])
        seg_boundary = sorted(seg_boundary)
        query = zip(seg_boundary[0: len(seg_boundary)], seg_boundary[1: len(seg_boundary)])

        for start, end in query:
            overlap = find_overlaps(start, end , tree)
            print '%s\t%s\t%s\t%s' % (chrom, start, end, len(overlap))     

if __name__ == '__main__':
    args = docopt.docopt(__doc__)
    chrom_lengths = find_chromosome_length(args['<seg_files>'])
    clogr_thresh = float(args['--clogr_thresh'])
    baf_thresh = float(args['--baf_thresh'])
    loh_regions_chromosome_wise = find_cnloh(args['<seg_files>'], clogr_thresh, baf_thresh)
    main(loh_regions_chromosome_wise, chrom_lengths)