#!/usr/bin/env python


import argparse
import mygene
import sys


"""
This script is meant to be run on the RefSeq GTF (limited format) downloaded from the UCSC table browser for hg38.
Using the mygene package, this script will query the gene_id listed for the gene name, gene description, and the numeric RefSeq geneid, and add these in as the last three columns.
"""


parser = argparse.ArgumentParser()
parser.add_argument('--gtf_in', required=True)
parser.add_argument('--gtf_out', required=True)
args = parser.parse_args()

fh = open(args.gtf_in, 'r')


#gene_id_list of uniq_gene_ids
gene_ids_annotated = {}

mg = mygene.MyGeneInfo()

#make dictionary of queried gene_id in order to get the name
i=0
for line in fh:
	line = line.strip()
	line = line.split('\t')
	line[8] = line[8].replace('"', '')
	line[8] = line[8].replace(';', '')
	eight_split = line[8].split()

	#index 1 is the gene id
	gene_id = eight_split[1]
	if gene_id not in gene_ids_annotated.keys(): #only want unique gene_ids here
		i += 1
		gene_ids_annotated[gene_id] = {}
		query_out = mg.query(gene_id, scopes='refseq')

		try:
			gene_ids_annotated[gene_id] = query_out['hits'][0]
		except:
			gene_ids_annotated[gene_id] = 'No hit'

	else:
		continue

fh.close()



## Now loop over the original gtf again, and write out a new gtf with 2 additional columns, gene_symbol and gene_name
fh = open(args.gtf_in, 'r')
out_fh = open(args.gtf_out, 'w')

for line in fh:
	line = line.strip()
	line = line.split('\t')
	original_8 = line[8]
	line[8] = line[8].replace('"', '')
	line[8] = line[8].replace(';', '')
	eight_split = line[8].split()

	#index 1 is the gene is
	gene_id = eight_split[1]

	line[8] = original_8

	#write out to file
	out_fh.write('\t'.join(line))
	if gene_ids_annotated[gene_id] == 'No hit':
		out_fh.write('\tNo_hit\tNo_hit\tNo_hit\n')
	else:
		out_fh.write('\t'+gene_ids_annotated[gene_id]['symbol']+'\t'+gene_ids_annotated[gene_id]['name']+'\t'+gene_ids_annotated[gene_id]['_id']+'\n')
fh.close()
out_fh.close()


