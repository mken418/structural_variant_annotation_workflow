#!/usr/bin/env python

import argparse
import subprocess
import os
import re
import sys
import collections

"""
	-- This script will makes a bed file containing all of the unique variants in a group of samples, and who has them, for downstream processing.
	-- Input is a directory containing a VCF for each individual.
	-- The script won't properly work with VCFs containing multiple individuals.
	-- The script is not meant to work with breakpoint notation VCFs.
"""


parser = argparse.ArgumentParser()
parser.add_argument('--path', required=True, help="Path to dir containing VCFs for each individual in the group")
parser.add_argument('--min_size', type=int, default=0, help="minimum size of a variant to consider")
parser.add_argument('--pass_filter', default=False, action='store_true', help="pass in this flag to only include variants with PASS in the filter column of the VCF")
parser.add_argument('--output', required=True, help="prefix for output bed file")
args = parser.parse_args()




##################################################################################################################################
var_lens = {} #a dictionary containing all of the lengths of the variants
var_types = {} #a dictionary containing all of the types of variants
who_has_it = {} #a dictionary containing all of the individuals who have the variant

#parse over all the sample vcfs in the directory
for vcf in os.listdir(args.path):
	match = re.match(r'(.+)\.vcf', vcf)
	sample = match.group(1)

	fh = open(args.path+'/'+vcf, 'r')
	for line in fh:
		if line.startswith('#'):
			continue

		line = line.strip()
		line = line.split('\t')

		if args.pass_filter and line[6] != "PASS":
			continue

		size = abs(len(line[3]) - len(line[4]))
		if size < args.min_size:
			continue


		contig = line[0].replace('_', '-') #to prevent parsing mistakes later
		locus = contig + '_' + line[1] + '_' + str(int(line[1]) + (len(line[3])-1)) #locus is the ref allele size -1 to account for the first base location already in VCF

		if locus not in var_lens.keys():
			var_lens[locus] = []
		if size not in var_lens[locus]:
			var_lens[locus].append(size)

		#get rough var_type based on size difference of ref and alt allele
		if (len(line[3]) - len(line[4])) == 0:
			var_type = 'same_length_var'

		elif (len(line[3]) - len(line[4])) < 0:
			var_type = 'ins'

		elif (len(line[3]) - len(line[4])) > 0:
			var_type = 'del'

		if locus not in var_types.keys():
			var_types[locus] = []
		if var_type not in var_types[locus]:
			var_types[locus].append(var_type)

		if locus not in who_has_it.keys():
			who_has_it[locus] = []
		who_has_it[locus].append(sample)	
	fh.close()




#write out combined unique variant to file
out_fh = open(args.output+'.tsv', 'w')
out_fh.write('#contig\tstart_locus\tstop_locus\tlength\tvar_type\tindividuals\ttotal_individuals\n')

for key, value in who_has_it.items():
	match = re.match(r'(.+)_(.+)_(.+)', key)
	contig, start, stop = match.group(1), match.group(2), match.group(3)

	out_fh.write(contig+'\t'+start+'\t'+stop+'\t')

	#could be multiple lengths of alt allelles, considering as the same record since the ref size is the same though, will report all lengths found in that 
	#range though. If there is multiple lengths for different individuals at the locus, its not clear in the output file. The user can look back at the VCFs.
	lens = [str(x) for x in var_lens[key]]
	out_fh.write(','.join(lens)+'\t')

	#write out the var types present at that location. same logic as lengths
	out_fh.write(','.join(var_types[key])+'\t')

	out_fh.write(','.join(value)+'\t'+str(len(value))+'\n')
out_fh.close()


#sort the output file
proc = subprocess.Popen(['sort', '-k', '1,1', '-k', '2,2n'], stdin=open(args.output+'.tsv', 'r'), stdout=subprocess.PIPE)
proc_string = str(proc.communicate()[0], 'utf-8')
out_fh = open(args.output+'.sorted.tsv', 'w')
out_fh.write(proc_string)
out_fh.close()
