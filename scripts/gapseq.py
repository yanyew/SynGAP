#!/usr/bin/env python3
# -*- coding:utf-8 -*-
import sys
import os
import re
from Bio import SeqIO


def rc(str):
    seq = str.upper()
    rev_seq = seq[::-1]
    rc_seq = ''
    bp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    standard = ['A', 'T', 'C', 'G']
    for nt in rev_seq:
        if nt in standard:
            rc_seq += bp[nt]
        else:
            rc_seq += nt
    return rc_seq


def dict_fa(infile):
    inf = open(infile, 'r')
    seqs = {}
    seq_id = ''
    while True:
        line = inf.readline()
        if not line: break
        if line[0] == '>':
            seq_id = line[1:].strip()
            seqs[seq_id] = ''
        else:
            seqs[seq_id] += line.strip()
    inf.close()
    return seqs


def split_str(str, length):
    str_arr = re.findall(r'.{%d}' % int(length), str)
    return str_arr


def list_gapbed(infile):
    inf_gapbed = open(infile, 'r')
    gapbed_list = inf_gapbed.readlines()
    return gapbed_list


def dict_gapbed(gapbed):
    gap_list = gapbed.split(';')
    gap_type = {}
    gap_chr = {}
    gap_start = {}
    gap_end = {}
    for gap in gap_list:
        if gap != '':
            gap = gap.rstrip('\n')
            type = gap.split('\t')[0]
            chr = gap.split('\t')[1]
            start = gap.split('\t')[2]
            end = gap.split('\t')[3]
            gap_id = gap.split('\t')[4]
            gap_type[gap_id] = type
            gap_chr[gap_id] = chr
            gap_start[gap_id] = start
            gap_end[gap_id] = end
    return gap_type, gap_chr, gap_start, gap_end


def gapseq(ingenome_sp1, ingenome_sp2, inpep_sp1, inpep_sp2, ingapbed_sp1, ingapbed_sp2):
    ingenomef_sp1 = open(ingenome_sp1, 'r')
    ingenomef_sp2 = open(ingenome_sp2, 'r')
    print('[\033[0;36mINFO\033[0m] Load file `'  + ingenome_sp1 + '`')
    genomeseqs_sp1 = SeqIO.to_dict(SeqIO.parse(ingenomef_sp1, 'fasta'))
    print('[\033[0;36mINFO\033[0m] Load file `'  + ingenome_sp2 + '`')
    genomeseqs_sp2 = SeqIO.to_dict(SeqIO.parse(ingenomef_sp2, 'fasta'))
    print('[\033[0;36mINFO\033[0m] Load file `'  + inpep_sp1 + '`')
    pepseqs_sp1 = dict_fa(inpep_sp1)
    print('[\033[0;36mINFO\033[0m] Load file `'  + inpep_sp2 + '`')
    pepseqs_sp2 = dict_fa(inpep_sp2)
    print('[\033[0;36mINFO\033[0m] Load file `'  + ingapbed_sp1 + '`')
    gapbedsp1_list = list_gapbed(ingapbed_sp1)
    print('[\033[0;36mINFO\033[0m] Load file `'  + ingapbed_sp2 + '`')
    gapbedsp2_list = list_gapbed(ingapbed_sp2)
    os.mkdir('./query')
    os.mkdir('./subject')
    outf_querynumber = open('./query/querynumber.txt', 'w')
    outf_subjectnumber = open('./subject/subjectnumber.txt', 'w')
    outf_querynumber.write(str(len(gapbedsp1_list)))
    outf_subjectnumber.write(str(len(gapbedsp2_list)))
    print('[\033[0;36mINFO\033[0m] Preparing ...')
    for n in range(len(gapbedsp1_list)):
        gaptype_sp1, gapchr_sp1, gapstart_sp1, gapend_sp1 = dict_gapbed(gapbedsp1_list[n])
        gaptype_sp2, gapchr_sp2, gapstart_sp2, gapend_sp2 = dict_gapbed(gapbedsp2_list[n])
        outf_query = open('./query/' + str(n) + '.query.fa', 'a')
        outf_subject = open('./subject/' + str(n) + '.subject.fa', 'a')
        for gap_id_sp1 in gaptype_sp1:
            if gaptype_sp1[gap_id_sp1] == '1':
                if gap_id_sp1 in pepseqs_sp1:
                    gapseq_sp1 = pepseqs_sp1[gap_id_sp1]
                    outf_query.write('>' + gap_id_sp1 + '\n' + str(gapseq_sp1) + '\n')
            elif gaptype_sp1[gap_id_sp1] == '2':
                chr_sp1 = genomeseqs_sp1[gapchr_sp1[gap_id_sp1]]
                gapseq_sp1 = chr_sp1.seq[int(gapstart_sp1[gap_id_sp1]):int(gapend_sp1[gap_id_sp1])]
                outf_subject.write('>' + gap_id_sp1 + '\n' + str(gapseq_sp1) + '\n')
        for gap_id_sp2 in gaptype_sp2:
            if gaptype_sp2[gap_id_sp2] == '1':
                chr_sp2 = genomeseqs_sp2[gapchr_sp2[gap_id_sp2]]
                gapseq_sp2 = chr_sp2.seq[int(gapstart_sp2[gap_id_sp2]):int(gapend_sp2[gap_id_sp2])]
                outf_subject.write('>' + gap_id_sp2 + '\n' + str(gapseq_sp2) + '\n')
            elif gaptype_sp2[gap_id_sp2] == '2':
                if gap_id_sp2 in pepseqs_sp2:
                    gapseq_sp2 = pepseqs_sp2[gap_id_sp2]
                    outf_query.write('>' + gap_id_sp2 + '\n' + str(gapseq_sp2) + '\n')
        outf_query.close()
        outf_subject.close()
        outf_querynumber.close()
        outf_subjectnumber.close()


if __name__ == '__main__':
    ingenome_sp1 = sys.argv[1]
    ingenome_sp2 = sys.argv[2]
    inpep_sp1 = sys.argv[3]
    inpep_sp2 = sys.argv[4]
    ingapbed_sp1 = sys.argv[5]
    ingapbed_sp2 = sys.argv[6]
    gapseq(ingenome_sp1, ingenome_sp2, inpep_sp1, inpep_sp2, ingapbed_sp1, ingapbed_sp2)