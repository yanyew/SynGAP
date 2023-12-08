#!/usr/bin/env python3
# -*- coding:utf-8 -*-
import sys


def dict_modifiedgff(infile):
    inf_modifiedgff = open(infile, 'r')
    transcripts = []
    mRNAs = {}
    miss_annotated_mRNAs = []
    while True:
        line = inf_modifiedgff.readline()
        if not line: break
        type = line.split('\t')[2]
        attributes = line.split('\t')[8].rstrip('\n')
        attributes_list = attributes.split(';')
        attributes_dict = {}
        for a in attributes_list:
            if a != '':
                attributes_dict[a.split('=')[0]] = a.split('=')[1]
        if type == 'mRNA':
            name = attributes_dict['ID']
            transcripts.append(name)
            if name not in mRNAs:
                mRNAs[name] = ''
                mRNAs[name] += line
            elif name in mRNAs:
                mRNAs[name] += line
            if attributes_dict['polish_type'] == 'miss_annotated':
                miss_annotated_mRNAs.append(name)
        elif type == 'CDS':
            name = attributes_dict['Parent']
            if name in mRNAs:
                mRNAs[name] += line
            elif name not in mRNAs:
                mRNAs[name] = ''
                mRNAs[name] += line
    inf_modifiedgff.close()
    return transcripts, mRNAs, miss_annotated_mRNAs


def polishtypescan(modified_gff, miss_annotated_gff, mis_annotated_gff):
    print('[\033[0;36mINFO\033[0m] Load file `\033[0;35m' + modified_gff + '\033[0m`')
    transcripts, mRNAs, miss_annotated_mRNAs = dict_modifiedgff(modified_gff)
    outf_miss_annotated_gff = open(miss_annotated_gff, 'w')
    outf_mis_annotated_gff = open(mis_annotated_gff, 'w')
    print('[\033[0;36mINFO\033[0m] Scanning ...')
    for name in mRNAs:
        if name in miss_annotated_mRNAs:
            outf_miss_annotated_gff.write(mRNAs[name])
        elif name not in miss_annotated_mRNAs:
            outf_mis_annotated_gff.write(mRNAs[name])
    outf_miss_annotated_gff.close()
    outf_mis_annotated_gff.close()


if __name__ == '__main__':
    modified_gff = sys.argv[1]
    miss_annotated_gff = sys.argv[2]
    mis_annotated_gff = sys.argv[3]
    polishtypescan(modified_gff, miss_annotated_gff, mis_annotated_gff)
