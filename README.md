# structural_variant_annotation_workflow
This set of scripts can be used take individual VCFs from one or more groups of samples, and determine the number of variants in each gene and gene region for each group. The final output file can be parsed for statistical analysis.

## Workflow order
1. aggregate_vcfs_into_group_bed.py
  * Run this script on an input directory of individual VCFs to combine all the variants into a bed-style TSV file containing each unique variant, variant type, and the individuals who have it.
  * Input VCFs must be individual VCFs, and must be in the same directory. If you plan to run this script on multiple groups of samples, each group should be i na separate directory.

2. annotate_UCSC_RefSeq_gtf.py
  * This script will add 3 additional fields to the RefSeq GTF (limited format) file available for download at the UCSC table browser: https://genome.ucsc.edu/cgi-bin/hgTables.
  * It is necessary to download the GTF and annotate with this script before running step 3.
  
3. group_bed_by_genes.py
  * This script will take each group bed-style tsv file produced in step1, and determine the gene and gene regions each variant overlaps with from the annotated RefSeq GTF (limited format) produced in step 2.
  * This script will output a sorted file that shows the group couts for each gene region (introns, exons, and including 50kb upstream/downtream), and each whole gene window, and a separate tsv file for each group that maps the variant to the gene region it was in.
