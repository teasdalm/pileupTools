#! /bin/sh

# bim2table.sh
# Script to convert bim to interval_list file

# Filter file
awk 'FS = "\t" {print "chr"$1":"$4"-"$4"}' $1

