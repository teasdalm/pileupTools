#! /usr/bin/python

"""
pileupTools
A script to manipulate GATK bam files
Matthew Teasdale 2016
MIT licence
"""

import argparse
from collections import Counter


def find_max_base(alleles):
    return Counter(allele)[0][0]


def filter_line(pileup_file_list):
    # filter SNP name
    long_snp_name = pileup_file_list[6].split('\t')[2]
    short_snp_name = long_snp_name.replace(']', '')
    current_line = {'chrom': pileup_file_list[0],
                    'pos': pileup_file_list[1],
                    'ref_base': pileup_file_list[2],
                    'alleles': pileup_file_list[3],
                    'base_qualities': pileup_file_list[4],
                    'snp_name': short_snp_name
                    }
    return current_line


def parse_pileup_file(pileup_file):
    with open(pileup_file) as f:
        for each_record in f:
            pileup_cols = each_record.rstrip().split(" ")
            if pileup_ok(pileup_cols):
                line_info = filter_line(pileup_cols)
                print line_info


def pileup_ok(pileup_line_list):
    if len(pileup_line_list) == 7:
        return True
    elif len(pileup_line_list) == 6 and pileup_line_list[0] == '[REDUCE':
        print "NB: GATK pileup reporting reduced output only {} SNPs used".format(pileup_line_list[4])
        return False
    else:
        return False


def filter_bases():
    pass


def get_arguments():
    parser = argparse.ArgumentParser(description='Manipulate GATK pileup files')
    parser.add_argument(dest='filenames', metavar='filename', nargs='*')
    parser.add_argument('-q', metavar='Minimum base quality')

    return parser.parse_args()

if __name__ == '__main__':
    my_args = get_arguments()
    for each in my_args.filenames:
        parse_pileup_file(each)