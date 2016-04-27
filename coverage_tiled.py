#!/bin/env pypy
"""
A tool for extracting read_depth.
 
Usage:
  coverage_tiled tiled <bam_file> <window_size> [-s N] [-q N] [-S|-1|-2] [-P|-M] 

Options:
  -q --quality=N    Minimum alignment quality [default: 10].
  -s --step=N       Step size for window placement [default: window size / 2].
  -S --single       Use all reads for coverage calculation, not just paired.
  -P --plus         Calculate coverage only for the plus strand.
  -M --minus        Calculate coverage only for the minus strand.
  -h --help       Show this screen.
"""

from __future__ import print_function
import numpy as np
import docopt, sys, re
from pypette import shell_stdout
# from sam import read_sam, ref_sequence_sizes

def read_sam(sam_path, mode='', min_quality=0):
	view_options = ''
	flag_on = 0x0
	flag_off = 0x900       # Ignore secondary and supplementary alignments
	if 'a' in mode: flag_off |= 0x4                   # Aligned
	if 'A' in mode: flag_on |= 0x1; flag_off |= 0xc   # Both mates aligned
	if 'C' in mode: flag_on |= 0x3; flag_off |= 0xc   # Concordant read pairs
	if 'u' in mode: flag_on |= 0x4                    # Unaligned
	if '1' in mode: flag_on |= 0x40                   # First mates
	if '2' in mode: flag_on |= 0x80                   # Second mates
	if '+' in mode: flag_off |= 0x10                  # Plus strand only
	if '-' in mode: flag_on |= 0x10                   # Minus strand only
	if not 'D' in mode: flag_off |= 0x400             # Flagged duplicates
	
	view_options += '-f 0x%x -F 0x%x ' % (flag_on, flag_off)
	
	if min_quality > 0: view_options += '-q%d ' % min_quality
	out = shell_stdout('samtools view %s %s' % (view_options, sam_path))
	for line in out:
		yield line.decode('utf8').split('\t')

def ref_sequence_sizes(sam_path):
	out = shell_stdout('samtools view -H %s' % sam_path)
	chr_sizes = {}
	for line in out:
		m = re.match('@SQ\tSN:(\w+)\tLN:(\d+)', line.decode('utf8'))
		if m:
			chr_sizes[m.group(1)] = int(m.group(2))
	return chr_sizes

def coverage_tiled(bam_path, window_size, quality, mode, step, max_frag_len=1000):
	#print(bam_path, window_size, quality, mode, step)
	#sys.exit()
	chr_sizes = ref_sequence_sizes(bam_path)
	chr_cov = { chr: np.zeros(int(size / step) + 1, np.uint32) for chr, size in chr_sizes.items() }
	
	win_overlap = window_size - step
	empty = np.array([])
	cov = empty
	chr = ''
	
	if 'A' in mode:
		# Count the entire length of the fragment.
		for al in read_sam(bam_path, mode, min_quality=quality):
			pos = int(al[3]); mpos = int(al[7])
			if al[6] != '=' or abs(pos - mpos) > max_frag_len: continue

			# FIXME: Spliced reads are tallied as if they were not spliced.
			start =     (min(pos, mpos) - win_overlap - 1) / step    
			stop =   (max(pos, mpos) + len(al[9]) - 1) / step    
			if al[2] != chr:
				chr = al[2]
				cov = chr_cov.get(al[2], empty)
			if not cov.size: continue
			start = max(start, 0)
			stop = min(stop, cov.size-1)
			cov[start:stop+1] += 1
		
	else:
		# Count all individual reads.
		for al in read_sam(bam_path, mode, min_quality=quality):
			pos = int(al[3])

			# FIXME: Spliced reads are tallied as if they were not spliced.
			start = (pos - win_overlap - 1) / step
			stop = (pos + len(al[9]) - 1) / step
			if al[2] != chr:
				chr = al[2]
				cov = chr_cov.get(al[2], empty)
			if not cov.size: continue
			start = max(start, 0)
			stop = min(stop, cov.size-1)
			cov[start:stop+1] += 1
	for chr in chr_cov:
		print('fixedStep chrom=%s start=%d step=%d' % (chr, step, step))
		for x in chr_cov[chr]: print(x)


if __name__ == '__main__':
	args = docopt.docopt(__doc__)
	if args['tiled']:
		wsize = int(args['<window_size>'])
		step = wsize / 2
		if args['--step'].isdigit(): step = int(args['--step'])
		
		mode = 'A'
		if '-S' in args: mode = 'a'
		if '-1' in args: mode = 'a1'
		if '-2' in args: mode = 'a2'
		if args['--plus']: mode += '+'
		if args['--minus']: mode += '-'
		
		coverage_tiled(args['<bam_file>'], wsize, quality=int(args['--quality']), mode=mode, step=step)