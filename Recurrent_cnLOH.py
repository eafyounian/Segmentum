#!/bin/env python
"""
A tool for extracting recurrent copy neutral LOH regions across one or many samples.

Usage:
  recurrent_cnLOH.py <seg_files>... [-c N] [-b N] 

Options:
  -h --help         Show this screen.
  -c --clogr_thresh=N   Coverage logratio must be below this threshold to call a copy neutral LOH region [default: 0.1]
  -b --baf_thresh=N     B-allele fraction must be below this threshold to call a copy neutral LOH region [default: 0.15]   
  
Author: Ebrahim Afyounian <ebrahim.afyounian@staff.uta.fi>
"""
from __future__ import print_function
import docopt, itertools, sys, math, random
from operator import itemgetter

class IntervalNode(object):
    def __init__(self, start, end, linenum=0, other=None):
        self.priority = math.ceil((-1.0 / math.log(.5)) * math.log( -1.0 / (random.uniform(0,1) - 1)))
        self.start = start
        self.end = end
        self.maxend = self.end
        self.minend = self.end
        self.left = None
        self.right = None
        self.linenum = linenum
        self.other = other
    def insert(self, start, end, linenum=0, other=None):
        root = self
        if start > self.start:
            # insert to right tree
            if self.right:
                self.right = self.right.insert(start, end, linenum, other)
            else:
                self.right = IntervalNode(start, end, linenum, other)
            # rebalance tree
            if self.priority < self.right.priority:
                root = self.rotateleft()
        else:
            # insert to left tree
            if self.left:
                self.left = self.left.insert(start, end, linenum, other)
            else:
                self.left = IntervalNode(start, end, linenum, other)
            # rebalance tree
            if self.priority < self.left.priority:
                root = self.rotateright()
        if root.right and root.left: 
            root.maxend = max(root.end, root.right.maxend, root.left.maxend)
            root.minend = min(root.end, root.right.minend, root.left.minend)
        elif root.right: 
            root.maxend = max(root.end, root.right.maxend)
            root.minend = min(root.end, root.right.minend)
        elif root.left:
            root.maxend = max(root.end, root.left.maxend)
            root.minend = min(root.end, root.left.minend)
        return root

    def rotateleft( self ):
        root = self.right
        self.right = self.right.left
        root.left = self
        if self.right and self.left: 
            self.maxend = max(self.end, self.right.maxend, self.left.maxend)
            self.minend = min(self.end, self.right.minend, self.left.minend )
        elif self.right:
            self.maxend = max(self.end, self.right.maxend)
            self.minend = min(self.end, self.right.minend)
        elif self.left:
            self.maxend = max(self.end, self.left.maxend)
            self.minend = min(self.end, self.left.minend )
        return root

    def rotateright( self ):
        root = self.left
        self.left = self.left.right
        root.right = self
        if self.right and self.left: 
            self.maxend = max(self.end, self.right.maxend, self.left.maxend)
            self.minend = min(self.end, self.right.minend, self.left.minend )
        elif self.right:
            self.maxend = max(self.end, self.right.maxend)
            self.minend = min(self.end, self.right.minend)
        elif self.left:
            self.maxend = max(self.end, self.left.maxend)
            self.minend = min(self.end, self.left.minend )
        return root

    def intersect(self, start, end, report_func ):
        if start < self.end and end > self.start: report_func( self)
        if self.left and start < self.left.maxend:
            self.left.intersect(start, end, report_func)
        if self.right and end > self.start:
            self.right.intersect(start, end, report_func)

def find_chromosome_length(seg_files):
    """finds the chromosome length from the .seg file and returns a dictionary of lengths"""
    chrom_lengths = {}
    chrom = {}
    
    seg_file = open(seg_files[0], 'r')
    # seg_file = open(seg_files[2], 'r')          ## just a Quick fix for a time when one chrom is missing. This should be fixed properly. 13.01.2016
    header = seg_file.readline()

    for line in seg_file:
        key = line.strip().split('\t')[1]
        seg_end = int(line.strip().split('\t')[3])
        if key not in chrom:
            chrom[key] = [seg_end]
        elif key in chrom:
            chrom[key].append(seg_end)
    seg_file.close()
    
    for key in chrom:
        chrom[key] = sorted(chrom[key])
        chrom_lengths[key] = chrom[key][-1] 
    
    return chrom_lengths

def find_overlaps(start, end, tree):
    """returns overlapping intervals from a interval tree"""
    out = []
    tree.intersect( start, end, lambda x: out.append(x) )
    return [(x.start, x.end) for x in out]
    
def find_cnloh(seg_files, clogr_thresh, baf_thresh):
    """finds regions with loss of heterozygozity in from several .seg files and returns the aggregated loh regions"""
    loh_regions_chromosome_wise = {}
    loh_regions = []
    
    for seg_file in seg_files:
        file = open(seg_file, 'r')
        header = file.readline().split('\t')
        for line in file:
            cols = line.strip().split('\t')
            if abs(float(cols[4])) <= clogr_thresh and float(cols[5]) >= (0.5 - baf_thresh):
                loh_regions.append(cols[0:4])
        file.close()
    
    loh_regions = sorted(loh_regions, key=itemgetter(1)) ##sorting based on chromosome name
    
    for region in loh_regions:
        if region[1] not in loh_regions_chromosome_wise:
            loh_regions_chromosome_wise[region[1]] = [(int(region[2]),int(region[3]))]
        elif region[1] in loh_regions_chromosome_wise:
            loh_regions_chromosome_wise[region[1]].append((int(region[2]),int(region[3])))

    for chrom in loh_regions_chromosome_wise:
        loh_regions_chromosome_wise[chrom] = sorted(loh_regions_chromosome_wise[chrom], key=itemgetter(0))
    return loh_regions_chromosome_wise

def find_recurrent(loh_regions_chromosome_wise, chrom_lengths):

    print('#chrom\tseg.start\tseg.end\toccurrence_number')
    for chrom in loh_regions_chromosome_wise:
        regions = loh_regions_chromosome_wise[chrom]

        ## creating the interval tree
        start, end = regions[0]
        tree = IntervalNode(start, end)
        for start, end in regions[1:]:
            tree = tree.insert(start, end)

        ## generating the query 
        seg_boundary =  set(itertools.chain(*regions))
        seg_boundary.add(0)
        seg_boundary.add(chrom_lengths[chrom])
        seg_boundary = sorted(seg_boundary)
        query = zip(seg_boundary[0: len(seg_boundary)], seg_boundary[1: len(seg_boundary)])

        for start, end in query:
            overlap = find_overlaps(start, end , tree)
            if len(overlap) == 0: continue
            print('%s\t%s\t%s\t%s' % (chrom, start, end, len(overlap)))

if __name__ == '__main__':
    args = docopt.docopt(__doc__)
    chrom_lengths = find_chromosome_length(args['<seg_files>'])
    clogr_thresh = float(args['--clogr_thresh'])
    baf_thresh = float(args['--baf_thresh'])
    loh_regions_chromosome_wise = find_cnloh(args['<seg_files>'], clogr_thresh, baf_thresh)
    find_recurrent(loh_regions_chromosome_wise, chrom_lengths)