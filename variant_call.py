#!/bin/env pypy

"""
Tools for calling short nucleotide variants.

Usage:
  variant_call calculate BAF <genome_fasta> <SNP_position> <tumor> <normal> [--hetz=N:R] [-q N] [-r REGION]

Options:
  -h --help         Show this screen
  -r <region>       Restrict analysis to chromosomal region
  -q --quality=N    Minimum mapping quality score [default: 10]
  --hetz=N:R        Minimum evidence for heterozygous [default: 4:0.25]
"""

from __future__ import print_function
import sys, docopt, re, os
from collections import defaultdict
from pypette import zopen, shell, shell_stdout, info, error

gt_symbols = ['', '0/0', '0/1', '1/1']

def simple_pileup(bam_paths, genome_path, kgenomes_path, min_mapq=10, min_alt_alleles=3,
	region=None):
	
	helper_dir = os.path.dirname(os.path.realpath(__file__)) + '/compiled'
	
	options = []
	if region:
		options.append('%s %s' % ('-l' if region.endswith('.bed') else '-r', region))
	
	# samtools mpileup will automatically ignore alignments flagged as
	# duplicates
	cmd = 'samtools mpileup -d 100000 -A -x -R -sB %s -q0 -l %s -f %s %s | %s/spileup %d %d' % (' '.join(options), kgenomes_path, genome_path,
		' '.join(bam_paths), helper_dir, min_alt_alleles, min_mapq)
	#info('Pre-filtering mutations with the following command:\n%s' % cmd)
	return shell_stdout(cmd)


def call_genotypes(alt, total, options):
	# 0 = unknown, 1 = ref, 2 = hetz, 3 = homz
	gtypes = [0] * len(alt)
	for s in range(len(alt)):
		if total[s] == 0: continue
		if total[s] - alt[s] >= options.min_ref_reads and \
			(total[s] - alt[s]) / total[s] >= options.min_ref_ratio:
			gtypes[s] = 1
		ratio = float(alt[s]) / total[s]
		if alt[s] >= options.min_hetz_reads and ratio >=options.min_hetz_ratio:
			gtypes[s] = 2
		if alt[s] >= options.min_homz_reads and ratio >=options.min_homz_ratio:
			gtypes[s] = 3
	return gtypes


################
# VARIANT CALL #
################

def calculate_BAF(bam_paths, genome_path, kgenomes_path, options):
	#print(bam_paths, genome_path, options.region, options.homz)
	gt_symbols = ['', '0/0', '0/1', '1/1']
	if not os.path.exists(genome_path):
		error('Could not find genome FASTA file %s.' % genome_path)

	if options.region:
		for bam_path in bam_paths:
			if not os.path.exists(bam_path + '.bai'):
				error('No index found for BAM file %s.' % bam_path)
	
	samples = [os.path.basename(p).replace('.bam', '') for p in bam_paths]

	# print('CHROM\tPOSITION\tREF\tALT\t%s' % '\t'.join(samples))
	print('CHROM\tPOSITION\tREF\tALT\t%s' %samples[0])

	ignore_mapq = [False] * len(samples)
	if options.ignore_mapq:
		for s, sample in enumerate(samples):
			if re.search(options.ignore_mapq, sample) != None:
				ignore_mapq[s] = True
				info('Ignoring mapping quality for sample %s.' % sample)
	
	for line in simple_pileup(bam_paths, genome_path, kgenomes_path,
		min_mapq=options.min_mapq, min_alt_alleles=(0 if options.keep_all else options.min_hetz_reads),
		region=options.region):

		if type(line) == bytes:
			line = line.decode('utf8')

		tokens = line[:-1].split('\t')
		if len(tokens) < 3: error('Invalid spileup line:\n%s' % line)
		if tokens[2] == 'N': continue
		pileups = [p.split(' ') for p in tokens[3:]]
		#total_reads = np.zeros(len(samples))
		#allele_reads = defaultdict(lambda: np.zeros(len(samples)))

		total_reads = [0] * len(samples)
		allele_reads = defaultdict(lambda: [0] * len(samples))

		for s, pileup in enumerate(pileups):
			if len(pileup) < 3: continue
			for a in range(0, len(pileup), 3):
				count = int(pileup[a+1]) + \
					(int(pileup[a+2]) if ignore_mapq[s] else 0)
				total_reads[s] += count
				if pileup[a] != '.': allele_reads[pileup[a]][s] = count		

		# Call genotypes for each allele.
		# for alt, reads in allele_reads.iteritems():
		for alt, reads in allele_reads.items():
			genotypes = call_genotypes(reads, total_reads, options)

			# if not options.keep_all and all(gt < 2 for gt in genotypes): continue
			# if all(gt != 2 for gt in genotypes): continue
			if genotypes[1] != 2: continue
			
			gtypes = ('%s:%d:%d' % (gt_symbols[g], reads[s], total_reads[s])
				for s, g in enumerate(genotypes))
			# Reformat indels in VCF4 format
			ref = tokens[2]
			if len(alt) >= 2:
				if alt[1] == '+':    # Insertion
					alt = (ref if alt[0] == '.' else alt[0]) + alt[2:]
				elif alt[1] == '-':  # Deletion
					ref += alt[2:]
					alt = (ref[0] if alt[0] == '.' else alt[0])
			
			#######################
			## Hetrozygous bases ##
			#######################
			
			gt_list = list(gtypes)
			gt_col = gt_list[1] ## genotype for the normal sample
			genotype = gt_symbols.index(gt_col[:gt_col.find(':')])
			total_read = float(gt_col.split(':')[2])
			if not (genotype == 2 and total_read >= 15): continue
			
			#########################
			## calculating the BAF ##
			#########################
			
			read = gt_list[0].split(':')[1:3] ## reads for the tumor sample
			sys.stdout.write('\t'.join([tokens[0], tokens[1], ref, alt.upper()]))
			alt, total = float(read[0]), int(read[1])
			sys.stdout.write('\tNaN' if total == 0 else '\t%.2f' % (alt / total))
			sys.stdout.write('\n')


#######################
# COMMAND LINE PARSER #
#######################

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

if __name__ == '__main__':
	args = docopt.docopt(__doc__)
	if args['calculate'] and args['BAF']:
		options = Options(args)
		bam_files = [args['<tumor>'], args['<normal>']]
		calculate_BAF(bam_files, args['<genome_fasta>'], args['<SNP_position>'], options)
