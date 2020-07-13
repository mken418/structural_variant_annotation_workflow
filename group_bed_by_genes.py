#!/usr/bin/env python

import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument('--gtf', required=True, help='use the gff with the accession names converted to the chromosome names')
parser.add_argument('--input_bed', action='append', required=True, help='input bed files to find their location in the gtf (multiple entries are allowed)')
parser.add_argument('--sample_name', action='append', required=True, help='Sample names that are associated with each bed. These must be in the same order as the input bed files, and the same amount.')
args = parser.parse_args()


print(args)
print(len(args.input_bed), len(args.sample_name))

#####################################
#####Pre-flight checks###############
#####################################
if len(args.input_bed) != len(args.sample_name):
	print('There are not the same number of input bed files and sample names. Please ensure that there is a sample name passed in with the same order as the input bed file')
	sys.exit()


######################################
####GTF into dict parsing function####
######################################
def add_count_values_to_gene_dict(gene_dict, chromosome, gene_region, start, stop, strand, gene_name, gene_desc, transcript_id, uniq_gene_id, first_occurance_gene_id, first_occurance_transcript_id=True):
	"""
	This function WILL modify the dict passed in.
	Developed to reduce redundant code when building the gene dict
	first_occurance is a true/false boolean, necessary for how exons are counted
	since multiple transcripts belong to a gene_id, if first_occurance_gene_id == True, then first_occurance_transcript_id must == True
	setting first_occurance_transcript_id as a kwarg so I don't have to type in when first_occurance_gene_is is used, even though it goes unused in that case.
	"""
	if first_occurance_gene_id: #need to initialize things so far
		gene_dict[uniq_gene_id] = {}
		gene_dict[uniq_gene_id]['chromo'] = chromosome
		gene_dict[uniq_gene_id]['strand'] = strand

		gene_dict[uniq_gene_id]['transcripts'] = []

		gene_dict[uniq_gene_id]['name'] = gene_name
		gene_dict[uniq_gene_id]['product'] = gene_desc

		#make name list for each sample which will store the individuals who have the variant
		for sample in args.sample_name:
			gene_dict[uniq_gene_id][sample] = []
			#exon only counts
			gene_dict[uniq_gene_id]['{}_exons'.format(sample)] = []
		
		gene_dict[uniq_gene_id]['coords'] = []		
	

	if first_occurance_transcript_id:
		gene_dict[uniq_gene_id]['transcripts'].append(transcript_id)
		gene_dict[uniq_gene_id][transcript_id] = {}

		if gene_region == 'start_codon':
			gene_dict[uniq_gene_id][transcript_id]['start_codon'] = (start, stop)
			gene_dict[uniq_gene_id][transcript_id]['num_exons_so_far'] = 0

		#ignoring CDS for now, it represents non-UTR regions, im still interested in exons even if they are UTR
		elif gene_region == 'CDS':
			pass

		elif gene_region == 'stop_codon':
			gene_dict[uniq_gene_id][transcript_id]['stop_codon'] = (start, stop)
			gene_dict[uniq_gene_id][transcript_id]['num_exons_so_far'] = 0

		elif gene_region == 'exon':
			if gene_dict[uniq_gene_id]['strand'] == '+':
				gene_dict[uniq_gene_id][transcript_id]['exon_1'] = (start, stop)
				gene_dict[uniq_gene_id][transcript_id]['num_exons_so_far'] = 1
			elif gene_dict[uniq_gene_id]['strand'] == '-':
				gene_dict[uniq_gene_id][transcript_id]['exon_-1'] = (start, stop) #negative strand, (first exon occurance is technically the last exon)
				gene_dict[uniq_gene_id][transcript_id]['num_exons_so_far'] = 1
		else:
			print("Something;s wrong, the gene region is not one thats in a known GTF (because I checked!)")
			sys.exit()


	else:
		i = gene_dict[uniq_gene_id][transcript_id]['num_exons_so_far']
		if gene_region == 'start_codon':
			gene_dict[uniq_gene_id][transcript_id]['start_codon'] = (start, stop)
			gene_dict[uniq_gene_id][transcript_id]['num_exons_so_far'] = i

		elif gene_region == 'CDS':
			pass

		elif gene_region == 'stop_codon':
			gene_dict[uniq_gene_id][transcript_id]['stop_codon'] = (start, stop)
			gene_dict[uniq_gene_id][transcript_id]['num_exons_so_far'] = i

		elif gene_region == 'exon':
			if gene_dict[uniq_gene_id]['strand'] == '+':
				gene_dict[uniq_gene_id][transcript_id]['exon_{}'.format(i+1)] = (start, stop)
				gene_dict[uniq_gene_id][transcript_id]['num_exons_so_far'] = i+1
			elif gene_dict[uniq_gene_id]['strand'] == '-':
				gene_dict[uniq_gene_id][transcript_id]['exon_-{}'.format(i+1)] = (start, stop)
				gene_dict[uniq_gene_id][transcript_id]['num_exons_so_far'] = i+1

	gene_dict[uniq_gene_id]['coords'].append(int(start)) #the min and max of this list will be the start and end location of the gene, respectively
	gene_dict[uniq_gene_id]['coords'].append(int(stop))

	return


##########################################################################
#functions used when comparing variant loci to gene dict ("intersection")#
##########################################################################
def convert_neg_strand_num(num_exons, iteration, exon=True):
	if exon:
		return (-iteration % num_exons) + 1
	else: #its an intron
		return (-iteration % num_exons)


def check_in(start1, stop1, start2, stop2):
	if (start1 >= start2 and stop1 <= stop2) or (start1 <= start2 and stop1 >= start2) or (start1 <= stop2 and stop1 >= stop2):
		return True
	else:
		return False


def add_to_group_dict(gene_dict, group_dict, gene_id, transcript, gene_region, people, out_fh, line, gene_region_output_formatted):
	"""
	Called in compare_variant_and_gene_copy to reduce redundant code.
	Will modify group_dict!!
	called if the gene region has a group variant overlapping, adds to the group dict, and then writes out the region and original line from the variants bed file along with gene region
	"""
	if transcript not in group_dict[gene_id].keys():
		group_dict[gene_id][transcript] = {}

	if gene_region not in group_dict[gene_id][transcript].keys():
		group_dict[gene_id][transcript][gene_region] = {}
		group_dict[gene_id][transcript][gene_region]['people'] = []
		group_dict[gene_id][transcript][gene_region]['counts'] = 0

	j = 0
	for person in people:
		if person not in group_dict[gene_id][transcript][gene_region]['people']:
			group_dict[gene_id][transcript][gene_region]['people'].append(person)
			group_dict[gene_id][transcript][gene_region]['counts'] += 1
			if j == 0:
				out_fh.write('\t'.join(line)+'\t'+gene_dict[gene_id]['name']+'\t'+gene_dict[gene_id]['product']+'\t'+transcript+'\t'+gene_region_output_formatted+'\n')
		j += 1
	return


def compare_variant_and_gene(gene_dict, group_dict, gene_id, people, group, start, stop, out_fh, line):
	"""
	Called after checking whether group variants were within the gene region
	This checks each possible gene region and calls add_to_group_dict to modify group_dict if there are variants in those gene regions
	Will modify group dict though calls to add_to_group_dict
	"""
	if gene_id not in group_dict.keys():
		group_dict[gene_id] = {}
		group_dict[gene_id]['total'] = 0
		group_dict[gene_id]['total_exons'] = 0

	for person in people:
		if person not in gene_dict[gene_id][group]:
			gene_dict[gene_id][group].append(person)
			group_dict[gene_id]['total'] += 1

	for transcript in gene_dict[gene_id]['transcripts']:
		#first check upstream for this transcript
		if check_in(start, stop, int(gene_dict[gene_id][transcript]['upstream'][0]), int(gene_dict[gene_id][transcript]['upstream'][1])):
			add_to_group_dict(gene_dict, group_dict, gene_id, transcript, 'upstream', people, out_fh, line, 'upstream')

		if check_in(start, stop, int(gene_dict[gene_id][transcript]['downstream'][0]), int(gene_dict[gene_id][transcript]['downstream'][1])):
			add_to_group_dict(gene_dict, group_dict, gene_id, transcript, 'downstream', people, out_fh, line, 'downstream')

		if gene_dict[gene_id]['strand'] == '+':
			i = 1
			while i <= gene_dict[gene_id][transcript]['num_exons_so_far']:
				if check_in(start, stop, int(gene_dict[gene_id][transcript]['exon_{}'.format(i)][0]), int(gene_dict[gene_id][transcript]['exon_{}'.format(i)][1])):
					add_to_group_dict(gene_dict, group_dict, gene_id, transcript, 'exon_{}'.format(i), people, out_fh, line, 'exon_{}'.format(i))

					#TO ONLY COUNT EXONS
					for person in people:
						if person not in gene_dict[gene_id]['{}_exons'.format(group)]:
							gene_dict[gene_id]['{}_exons'.format(group)].append(person)
							group_dict[gene_id]['total_exons'] += 1
						else:
							continue

				if i < gene_dict[gene_id][transcript]['num_exons_so_far']: #the number of introns will always be 1 less than this  
					if check_in(start, stop, int(gene_dict[gene_id][transcript]['intron_{}'.format(i)][0]), int(gene_dict[gene_id][transcript]['intron_{}'.format(i)][1])):
						add_to_group_dict(gene_dict, group_dict, gene_id, transcript, 'intron_{}'.format(i), people, out_fh, line, 'intron_{}'.format(i))

				i += 1

		elif gene_dict[gene_id]['strand'] == '-':
			i = 1
			while i <= gene_dict[gene_id][transcript]['num_exons_so_far']:
				if check_in(start, stop, int(gene_dict[gene_id][transcript]['exon_-{}'.format(i)][0]), int(gene_dict[gene_id][transcript]['exon_-{}'.format(i)][1])):
					add_to_group_dict(gene_dict, group_dict, gene_id, transcript, 'exon_-{}'.format(i), people, out_fh, line, 'exon_{}'.format(convert_neg_strand_num(gene_dict[gene_id][transcript]['num_exons_so_far'], i, True)))

					#TO ONLY COUNT EXONS
					for person in people:
						if person not in gene_dict[gene_id]['{}_exons'.format(group)]:
							gene_dict[gene_id]['{}_exons'.format(group)].append(person)
							group_dict[gene_id]['total_exons'] += 1
						else:
							continue

				if i < gene_dict[gene_id][transcript]['num_exons_so_far']: #the number of introns will always be 1 less than this
					if check_in(start, stop, int(gene_dict[gene_id][transcript]['intron_-{}'.format(i)][0]), int(gene_dict[gene_id][transcript]['intron_-{}'.format(i)][1])):
						add_to_group_dict(gene_dict, group_dict, gene_id, transcript, 'intron_-{}'.format(i), people, out_fh, line, 'intron_{}'.format(convert_neg_strand_num(gene_dict[gene_id][transcript]['num_exons_so_far'], i, False)))

				i += 1


	return 


#################################################################################################
#########Functions used when comparing the group dicts at the end to write out counts############
#################################################################################################
