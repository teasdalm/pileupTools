#! /usr/bin/python

"""
pileupTools
A script to manipulate GATK bam files
Matthew Teasdale 2016
BSD 3-Clause Licensed
"""

import argparse
from collections import Counter
import random


def allele_check(snp_tuple,allele):
    reverse_trans = {'A': 'T',
                     'G': 'C',
                     'C': 'G',
                     'T': 'A'}
    if allele == snp_tuple[0] or allele == snp_tuple[1]:
        return True
    elif reverse_trans[allele] == snp_tuple[0] or reverse_trans[allele] == snp_tuple[1]:
        return True
    else:
        return False


def find_max_base(alleles):
    allele_counts = Counter(alleles).most_common()
    max_alleles = []
    max = allele_counts[0][1]
    for each_allele in allele_counts:
        if each_allele[1] == max:
            max_alleles.append(each_allele[0])
    return random.choice(max_alleles)


def check_bases(line_dict, base_qual_cutoff):
    current_quals = line_dict['base_qualities']
    current_alleles = line_dict['alleles']
    new_alleles = ""
    new_quals = ""

    for i in range(0, len(current_quals)):
        # filter for base quality
        base_qual = ord(current_quals[i]) - 33
        if base_qual >= base_qual_cutoff:
            # Filter for none ATGC
            if current_alleles[i] in "AGTC":
                new_alleles = new_alleles + current_alleles[i]
                new_quals = new_quals + current_quals[i]
    line_dict['base_qualities'] = new_quals
    line_dict['alleles'] = new_alleles
    return line_dict


def filter_line(pileup_line_list):
    # filter SNP name
    info_cols = pileup_line_list[6].split('\t')
    short_snp_name = info_cols[2]
    snp_a = info_cols[3]
    snp_b = info_cols[4].rstrip(']')
    current_line = {'chrom': pileup_line_list[0],
                    'pos': pileup_line_list[1],
                    'ref_base': pileup_line_list[2],
                    'alleles': pileup_line_list[3].upper(),
                    'base_qualities': pileup_line_list[4],
                    'snp_name': short_snp_name,
                    'snps': (snp_a, snp_b)
                    }
    return current_line


def pileup_ok(pileup_line_list):
    if len(pileup_line_list) == 7:
        return True
    elif len(pileup_line_list) == 6 and pileup_line_list[0] == '[REDUCE':
        print("NB: GATK pileup reporting reduced output stating only {} SNPs used".format(pileup_line_list[5]))
        return False
    else:
        raise SystemExit('Bad pileup file!')


def parse_pileup_file(pileup_file, sample, base_min_quality):
    # dealing with files
    file_base = pileup_file.replace(".pileup", "")
    map_file = open(file_base + '.map', 'wt')
    ped_file = open(file_base + '.ped', 'wt')
    ped_file.write('{s} {s} 0 0 0 -9'.format(s=sample))

    # counter
    ct_good = 0
    ct_bad = 0

    # reading pileup
    with open(pileup_file) as f:
        for each_record in f:
            pileup_cols = each_record.rstrip().split(" ")
            # check pileup line
            if pileup_ok(pileup_cols):
                # filtering snps
                line_info = filter_line(pileup_cols)
                line_info = check_bases(line_info, base_min_quality)
                # check for missing data
                if len(line_info['alleles']) != 0:
                    ct_good += 1
                    allele_called = find_max_base(line_info['alleles'])
                    # check SNP matches ref alleles
                    if allele_check(line_info['snps'], allele_called):
                        # write to output files
                        ped_file.write(' {a} {a}'.format(a=allele_called))
                        map_file.write('{c}\t{s}\t0\t{p}\n'
                                       .format(c=line_info['chrom'],
                                               s=line_info['snp_name'],
                                               p=line_info['pos']))
                    else:
                        print('{} tri-allele'.format(line_info['snp_name']))
                        ct_bad += 1
                else:
                    ct_bad += 1
                    print('{} failed filters'.format(line_info['snp_name']))
    map_file.close()
    ped_file.write('\n')
    ped_file.close()

    print('Analysis finished wrote {:,} SNPs to output files removed {:,}.\nOutput PED ---> {}\nOutput MAP --->{}'
          .format(ct_good, ct_bad, file_base + '.ped', file_base + '.map'))


def get_arguments():
    parser = argparse.ArgumentParser(description='A tool to manipulate GATK pileup files')
    parser.add_argument(dest='filename', metavar='filename', help='Input pileup file [file.pileup].')
    parser.add_argument('-q',
                        metavar='Minimum Base Quality',
                        help='Minimum base quality to filter [30].',
                        dest='min_quality',
                        type=int,
                        default=30)
    parser.add_argument('-s',
                        metavar='Sample Name',
                        default='XXX',
                        type=str,
                        dest='sample_name',
                        help='Your sample name [XXX].')
    return parser.parse_args()


def main():
    my_args = get_arguments()
    parse_pileup_file(my_args.filename,
                      my_args.sample_name,
                      my_args.min_quality)
if __name__ == '__main__':
    main()
