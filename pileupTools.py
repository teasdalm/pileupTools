#! /usr/bin/python

"""
pileupTools
A script to manipulate GATK bam files
Matthew Teasdale 2016
MIT licence
"""

import argparse

def parse_pileup_file(pileup_file):
    pass

def filter_bases:
    pass

def get_arguments():
    parser = argparse.ArgumentParser(description='Manipulate GATK pileup files')
    parser.add_argument(dest='filenames', metavar='filename', nargs='*')
    parser.add_argument('-q', metavar='Minimum base quality')

    return parser.parse_args()

if __name__ == '__main__':
    my_args = get_arguments()
    for each in my_args.filenames:
        print each