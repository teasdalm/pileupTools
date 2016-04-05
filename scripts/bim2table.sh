#! /bin/sh

# bim2table.sh
# Script to convert bim to TABLE file 
# Print header
echo "HEADER	SNP-NAME	ALLELE-A	ALLELE-B"

# Filter file
awk 'FS = "\t" {print "chr"$1":"$4"-"$4"\t"$2"\t"$5"\t"$6}' $1 

