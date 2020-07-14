#!/usr/bin/env python

import argparse
import sys
import collections

parser = argparse.ArgumentParser()
parser.add_argument('--gtf', required=True, help='use the gff with the accession names converted to the chromosome names')
parser.add_argument('--input_bed', action='append', required=True, help='input bed files to find their location in the gtf (multiple entries are allowed)')
parser.add_argument('--sample_name', action='append', required=True, help='Sample names that are associated with each bed. These must be in the same order as the input bed files, and the same amount.')
parser.add_argument('--gene_group_output', required=True, help='prefix for the output file name, the script writes an unsorted output, then uses subprocess to output a sorted one')
args = parser.parse_args()


"""
Important notes:
	-- the GTF is in a custom format. So GTF based. This code is currently meant to parse the RefSeq GTF (limited format) downloaded from the UCSC table browser.
	This parsing code expects the GTF to have been queryed with the mygene python script to add the gene names, gene descriptions, and gene ids

	--the upstream and downstream regions for each gene are 50kb +/- the gene body

	--not including the alternate contigs for now. The random/Un/HLA contigs were parsed out when combining the variants for each group
"""


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
def check_group_prescence(gene_dict, gene_id, transcript, gene_region, gene_region_output_formatted, out_fh, sample_var_dicts):
	"""
	Checks if a variant present in one sample group is present in the other sample groups, and then writes out the group counts for each gene region into the output file
	"""
	out_fh.write(gene_dict[gene_id]['chromo']+'\t'+str(gene_dict[gene_id][transcript][gene_region][0])+'\t'+str(gene_dict[gene_id][transcript][gene_region][1])+'\t'+str(gene_id)+'\t'+gene_dict[gene_id]['strand']+'\t'+gene_dict[gene_id]['name']+'\t'+gene_dict[gene_id]['product']+'\t'+transcript+'\t'+gene_region_output_formatted+'\t')	

	group_counts = collections.OrderedDict()
	#key is sample group name, value is the group_dict
	for name,group_dict in sample_var_dicts.items():
		try: #will fail is gene region or transcript isnt present in the dict
			group_counts[name] = group_dict[gene_id][transcript][gene_region]['counts']
		except: 
			group_counts[name] = 0

	length = len(group_counts.keys())
	i = 0 
	for name, count in group_counts.items():
		i += 1
		if i < length:
			out_fh.write(str(count)+'\t')
		else:
			out_fh.write(str(count)+'\n')

	return	




##########################################################################################################################################
###########################################################    Code Body         #########################################################
##########################################################################################################################################

#####################
#parse gff into dict
gene_dict = {}

fh = open(args.gtf, 'r')
for line in fh:
	line = line.strip()
	line = line.split('\t')
	chromo = line[0]
	if 'alt' in chromo: #skipping alt contigs for now	
		continue

	gene_region = line[2]

	start = line[3]
	stop = line[4]

	strand = line[6]

	gene_name = line[9]
	gene_desc = line[10]

	if line[11] == "No_hit":
		continue

	uniq_gene_id = line[11]

	line[8] = line[8].replace('"', '')
	line[8] = line[8].replace(';', '')

	trans_col = line[8].split()
	transcript_id = trans_col[3]

	if uniq_gene_id not in gene_dict.keys():
		first_occurance_gene_id = True
		add_count_values_to_gene_dict(gene_dict, chromo, gene_region, start, stop, strand, gene_name, gene_desc, transcript_id, uniq_gene_id, first_occurance_gene_id)		

	elif chromo != gene_dict[uniq_gene_id]['chromo']: #uniq gene id is in the gene dictionary, but the chomosome positions do not match, happens on the X and Y chromosomes
		uniq_gene_id_chromo = uniq_gene_id + '_' + chromo #append to chromosome name to subsequent appearances of this gene_id

		if uniq_gene_id_chromo not in gene_dict.keys():
			first_occurance_gene_id = True
			add_count_values_to_gene_dict(gene_dict, chromo, gene_region, start, stop, strand, gene_name, gene_desc, transcript_id, uniq_gene_id_chromo, first_occurance_gene_id)

		else:
			first_occurance_gene_id = False
			if transcript_id not in gene_dict[uniq_gene_id_chromo].keys():
				first_occurance_transcript_id = True
				add_count_values_to_gene_dict(gene_dict, chromo, gene_region, start, stop, strand, gene_name, gene_desc, transcript_id, uniq_gene_id_chromo, first_occurance_gene_id, first_occurance_transcript_id)
			else: #transcript ID is in the unique gene IDs
				first_occurance_transcript_id = False
				add_count_values_to_gene_dict(gene_dict, chromo, gene_region, start, stop, strand, gene_name, gene_desc, transcript_id, uniq_gene_id_chromo, first_occurance_gene_id, first_occurance_transcript_id)

	else: #uniq gene id is in the dictionary so append the coords to main gene, and check the transcripts
		first_occurance_gene_id = False
		if transcript_id not in gene_dict[uniq_gene_id].keys():
			first_occurance_transcript_id = True
			add_count_values_to_gene_dict(gene_dict, chromo, gene_region, start, stop, strand, gene_name, gene_desc, transcript_id, uniq_gene_id, first_occurance_gene_id, first_occurance_transcript_id)
		else:
			first_occurance_transcript_id = False
			add_count_values_to_gene_dict(gene_dict, chromo, gene_region, start, stop, strand, gene_name, gene_desc, transcript_id, uniq_gene_id, first_occurance_gene_id, first_occurance_transcript_id)

fh.close()


################################################################################################################################################################################
#Now I need to loop over the newly created dictionary, to parse make the introns/exons, and upstream/downstream regions, that will then be checked along with parsing the genes
for gene_id in gene_dict.keys():
	for transcript in gene_dict[gene_id]['transcripts']:
		if gene_dict[gene_id]['strand'] == '+': #my introns/exons are going to be 1 indexed
			gene_dict[gene_id][transcript]['upstream'] = ( (int(gene_dict[gene_id][transcript]['exon_1'][0]) - 50000), (int(gene_dict[gene_id][transcript]['exon_1'][0]) - 1) ) #taking 50kb upstream as upstream
			last_exon_index = gene_dict[gene_id][transcript]['num_exons_so_far']
			gene_dict[gene_id][transcript]['downstream'] = ( (int(gene_dict[gene_id][transcript]['exon_{}'.format(last_exon_index)][1]) + 1), (int(gene_dict[gene_id][transcript]['exon_{}'.format(last_exon_index)][1]) + 50000) ) #taking 50kb downstream

			i = 1
			while i < gene_dict[gene_id][transcript]['num_exons_so_far']: #this number will actually be the total number of exons, since im looking for introns, will be 1 less
				#defining the gene regions here; regions are defined as tuple: (start, stop) introns between eoxns, etc
				#assuming genes start/end with exons, because that's how the gff represents them
				gene_dict[gene_id][transcript]['intron_{}'.format(i)] = ( (int(gene_dict[gene_id][transcript]['exon_{}'.format(i)][1]) + 1), (int(gene_dict[gene_id][transcript]['exon_{}'.format(i+1)][0]) -1 ) ) #Intron_i will be region from (exon_i + 1) to to beginning of (exon_i+1 -1)
				i += 1
		elif gene_dict[gene_id]['strand'] == '-': #have to think in terms of neg strand 3' to 5' reverse orientation
			first_exon_index = gene_dict[gene_id][transcript]['num_exons_so_far'] #will be the most negative number
			gene_dict[gene_id][transcript]['upstream'] = ( (int(gene_dict[gene_id][transcript]['exon_-{}'.format(first_exon_index)][1]) + 1), (int(gene_dict[gene_id][transcript]['exon_-{}'.format(first_exon_index)][1]) + 50000) )
			gene_dict[gene_id][transcript]['downstream'] = ( (int(gene_dict[gene_id][transcript]['exon_-1'][0]) - 50000), (int(gene_dict[gene_id][transcript]['exon_-1'][0]) - 1) )

			i=1
			while i < gene_dict[gene_id][transcript]['num_exons_so_far']: #this number will actually be the total number of exons, since im looking for introns, will be one less
				gene_dict[gene_id][transcript]['intron_-{}'.format(i)] = ( (int(gene_dict[gene_id][transcript]['exon_-{}'.format(i)][1]) + 1), (int(gene_dict[gene_id][transcript]['exon_-{}'.format(i+1)][0]) -1 ) )
				i += 1

	###now append to coords the upstream/downstream
	gene_dict[gene_id]['coords'].append(min(gene_dict[gene_id]['coords'])-50000)
	gene_dict[gene_id]['coords'].append(max(gene_dict[gene_id]['coords'])+50000)



########################################################
#Now parse the groups bed files against the dict
sample_var_dicts = collections.OrderedDict() #will be a dict of dictionaries
i = 0
for bed in args.input_bed:
	sample_var_dicts[args.sample_name[i]] = {}
	fh = open(bed, 'r')
	out_fh = open(args.args.sample_name[i]+'.tsv', 'w')

	for line in fh:
		if line.startswith('#'):
			continue

		line = line.strip()
		line = line.split('\t')

		chromo = line[0]
		#not including alt contigs for now
		if 'alt' in chromo:
			continue

		start = int(line[1])
		stop = int(line[2])

		people = line[5].split(',')	

		for gene_id in gene_dict.keys():
			if chromo == gene_dict[gene_id]['chromo']:
				if (start >= min(gene_dict[gene_id]['coords']) and stop <= max(gene_dict[gene_id]['coords'])) or (start <= min(gene_dict[gene_id]['coords']) and stop >= min(gene_dict[gene_id]['coords'])) or (start <= max(gene_dict[gene_id]['coords']) and stop >= max(gene_dict[gene_id]['coords'])): #if variant is  within the gene region including upstream
					compare_variant_and_gene(gene_dict, sample_var_dicts[args.sample_name[i]], gene_id, people, 'probands', start, stop, out_fh, line)
				else:
					continue
			else:
				continue

	fh.close()
	out_fh.close()



###########################################################
#write main output file with gene and gene region counts
out_fh = open(args.gene_group_output+'.tsv', 'w')
out_fh.write('#chromo\tstart\tstop\tgene_id\tstrand\tgene_name\tgene_product\ttranscript\tgene_region\t')
i = 0
for sample in args.sample_name:
	i += 1
	if i < len(args.sample_name):
		out_fh.write(sample+'\t')
	else:
		out_fh.write(sample+'\n')

#list will contain gene_ids and gene_id+'_'+transcript_id+'_'+region to account for regions that were added in previous group iterations
already_added = []

for sample in args.sample_name:
	for gene_id in sample_var_dicts[sample].keys():
		if gene_id not in already_added:
			out_fh.write(gene_dict[gene_id]['chromo']+'\t'+str(min(gene_dict[gene_id]['coords']))+'\t'+str(max(gene_dict[gene_id]['coords']))+'\t'+str(gene_id)+'\t'+gene_dict[gene_id]['strand']+'\t'+gene_dict[gene_id]['name']+'\t'+gene_dict[gene_id]['product']+'\t'+'full_window'+'\t'+'full_window'+'\t')
			i = 0
			for sample in args.sample_name:
				i += 1
				if i < len(args.sample_name):
					out_fh.write(str(len(gene_dict[gene_id][sample]))+'\t')
				else:
					out_fh.write(str(len(gene_dict[gene_id][sample]))+'\n')

			out_fh.write(gene_dict[gene_id]['chromo']+'\t'+str(min(gene_dict[gene_id]['coords']))+'\t'+str(max(gene_dict[gene_id]['coords']))+'\t'+str(gene_id)+'\t'+gene_dict[gene_id]['strand']+'\t'+gene_dict[gene_id]['name']+'\t'+gene_dict[gene_id]['product']+'\t'+'only_exon_counts'+'\t'+'only_exon_counts'+'\t')
			i = 0
			for sample in args.sample_name:
				i += 1
				if i < len(args.sample_name):
					out_fh.write(str(len(gene_dict[gene_id]['{}_exons'.format(sample)]))+'\t')
				else:
					out_fh.write(str(len(gene_dict[gene_id]['{}_exons'.format(sample)]))+'\n')

			already_added.append(gene_id)

		for transcript in gene_dict[gene_id]['transcripts']:
			if transcript not in sample_var_dicts[sample][gene_id].keys():
				continue

			if 'upstream' in sample_var_dicts[sample][gene_id][transcript].keys():
				if str(gene_id)+'_'+str(transcript)+'_'+'upstream' not in already_added:
					check_group_prescence(gene_dict, gene_id, transcript, 'upstream', 'upstream', out_fh, sample_var_dicts)
					already_added.append(str(gene_id)+'_'+str(transcript)+'_'+'upstream')

			if 'downstream' in sample_var_dicts[sample][gene_id][transcript].keys():
				if str(gene_id)+'_'+str(transcript)+'_'+'downstream' not in already_added:
					check_group_prescence(gene_dict, gene_id, transcript, 'downstream', 'downstream', out_fh, sample_var_dicts)
					already_added.append(str(gene_id)+'_'+str(transcript)+'_'+'downstream')

			if gene_dict[gene_id]['strand'] == '+':
				i = 0
				while i <= gene_dict[gene_id][transcript]['num_exons_so_far']:
					if 'exon_{}'.format(i) in sample_var_dicts[sample][gene_id][transcript].keys():
						if str(gene_id)+'_'+str(transcript)+'_'+'exon_{}'.format(i) not in already_added:
							check_group_prescence(gene_dict, gene_id, transcript, 'exon_{}'.format(i), 'exon_{}'.format(i), out_fh, sample_var_dicts)
							already_added.append(str(gene_id)+'_'+str(transcript)+'_'+'exon_{}'.format(i))

					if i < gene_dict[gene_id][transcript]['num_exons_so_far']:
						if 'intron_{}'.format(i) in sample_var_dicts[sample][gene_id][transcript].keys():
							if str(gene_id)+'_'+str(transcript)+'_'+'intron_{}'.format(i) not in already_added:
								check_group_prescence(gene_dict, gene_id, transcript, 'intron_{}'.format(i), 'intron_{}'.format(i), out_fh, sample_var_dicts)
								already_added.append(str(gene_id)+'_'+str(transcript)+'_'+'intron_{}'.format(i))

					i += 1

			elif gene_dict[gene_id]['strand'] == '-':
				i = 0
				while i <= gene_dict[gene_id][transcript]['num_exons_so_far']:
					if 'exon_-{}'.format(i) in sample_var_dicts[sample][gene_id][transcript].keys():
						if str(gene_id)+'_'+str(transcript)+'_'+'exon_-{}'.format(i) not in already_added:
							check_group_prescence(gene_dict, gene_id, transcript, 'exon_-{}'.format(i), 'exon_{}'.format(convert_neg_strand_num(gene_dict[gene_id][transcript]['num_exons_so_far'], i, True)), out_fh, sample_var_dicts)
							already_added.append(str(gene_id)+'_'+str(transcript)+'_'+'exon_-{}'.format(i))

					if i < gene_dict[gene_id][transcript]['num_exons_so_far']:
						if 'intron_-{}'.format(i) in sample_var_dicts[sample][gene_id][transcript].keys():
							if str(gene_id)+'_'+str(transcript)+'_'+'intron_-{}'.format(i) not in already_added:
								check_group_prescence(gene_dict, gene_id, transcript, 'intron_-{}'.format(i), 'intron_{}'.format(convert_neg_strand_num(gene_dict[gene_id][transcript]['num_exons_so_far'], i, False)), out_fh, sample_var_dicts)
								already_added.append(str(gene_id)+'_'+str(transcript)+'_'+'intron_-{}'.format(i))
					i += 1





out_fh.close()

#sort output to new file
proc = subprocess.Popen(['sort', '-k', '1,1', '-k', '2,2n'], stdin=open(args.gene_group_output+'.tsv', 'r'), stdout=subprocess.PIPE)
proc_string = str(proc.communicate()[0], 'utf-8')
out_fh = open(args.gene_group_output+'sorted.tsv', 'w')
out_fh.write(proc_string)
out_fh.close()
