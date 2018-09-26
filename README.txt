####################################################################################################################
####################################################################################################################
#### Bioinformatics technical challenge: Variant annotation
#### Author: Sunantha Sethuraman
#### Date: 09/07/2017
#### Use: Annotates the given vcf file with information from online resources
####################################################################################################################
####################################################################################################################

File provided: Challenge_data.vcf
Input file: Challenge_data_edited.vcf
Script: VariantAnnotation.R
Output file: AnnotatedTable.csv

Notes:
	1. The file provided was missing the "FILTER" information, which was necessary for bioconductor functions to read the VCF file. Hence the line "##FILTER=<ID=.,Description="Not Filtered">" was added to the provided file and the file was renamed as "Challenge_data_edited.vcf".
	2. No filters were applied during variant annotation.
	3. All annotations requested in the problem statement were performed.
	4. Annotation with the most deleterious possible outcome in case of multi-allelic positions was not possible due to a potential bug in the VariantFiltering package. In an ideal case, it would have involved annotation of an expandedVCF (in contrast to the collapsedVCF used in the code), followed by some simple filtering in the output dataframe.
