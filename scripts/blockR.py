#!/usr/bin/env python3
# -*- coding:utf-8 -*-
import os
import sys
import multiprocessing
import shlex
import subprocess
import numpy as np


def dict_bed(infile):
    inf_bed = open(infile, 'r')
    gene_bed = {}
    while True:
        line = inf_bed.readline()
        if not line: break
        chr = line.split('\t')[0]
        start = line.split('\t')[1]
        end = line.split('\t')[2]
        gene_id = line.split('\t')[3]
        if gene_id not in gene_bed:
            gene_bed[gene_id] = (chr, int(start), int(end))
    inf_bed.close()
    return gene_bed


def list_blocks(inf_anchors):
    in_anchors = open(inf_anchors, 'r')
    anchors = in_anchors.read()
    anchors_list = anchors.split('###\n')
    blocks_list = {}
    block_number = 0
    for block in anchors_list:
        if block != '':
            block_name = 'block' + str(block_number)
            blocks_list[block_name] = block.rstrip('\n')
            block_number += 1
        else:
            continue
    in_anchors.close()
    return blocks_list


def dict_block(block):
    block_line = {}
    line_num = 0
    lines = block.split('\n')
    while '' in lines:
        lines.remove('')
    for line in lines:
        sp1_geneid = line.split('\t')[0]
        sp2_geneid = line.split('\t')[1]
        block_line[line_num] = (sp1_geneid, sp2_geneid)
        line_num += 1
    return block_line


def dict_anchor(inf_anchors, infasta1, infasta2, inbed_sp1, inbed_sp2):
    blocks_list = list_blocks(inf_anchors)
    gene_bed_sp1 = dict_bed(inbed_sp1)
    gene_bed_sp2 = dict_bed(inbed_sp2)
    synteny_pair = {}
    num = 0
    block_range_sp1 = {}
    block_range_sp2 = {}
    block_start_sp1 = {}
    block_end_sp1 = {}
    block_start_sp2 = {}
    block_end_sp2 = {}
    for block_name in blocks_list:
        block_start_sp1[block_name] = []
        block_end_sp1[block_name] = []
        block_start_sp2[block_name] = []
        block_end_sp2[block_name] = []
        block = blocks_list[block_name]
        block_line = dict_block(block)
        for line in block_line:
            gene1 = block_line[line][0]
            gene2 = block_line[line][1]
            pair = (block_name, gene1, gene2)
            synteny_pair[num] = pair
            num += 1
            gene1_region = gene_bed_sp1[gene1]
            gene2_region = gene_bed_sp2[gene2]
            gene1_chr = gene1_region[0]
            gene1_start = gene1_region[1]
            gene1_end = gene1_region[2]
            gene2_chr = gene2_region[0]
            gene2_start = gene2_region[1]
            gene2_end = gene2_region[2]
            block_range_sp1[block_name] = [gene1_chr]
            block_range_sp2[block_name] = [gene2_chr]
            block_start_sp1[block_name].append(int(gene1_start))
            block_end_sp1[block_name].append(int(gene1_end))
            block_start_sp2[block_name].append(int(gene2_start))
            block_end_sp2[block_name].append(int(gene2_end))
        block_range_sp1[block_name].append(min(block_start_sp1[block_name]))
        block_range_sp1[block_name].append(max(block_end_sp1[block_name]))
        block_range_sp2[block_name].append(min(block_start_sp2[block_name]))
        block_range_sp2[block_name].append(max(block_end_sp2[block_name]))

    inf1 = open(infasta1, 'r')
    inf2 = open(infasta2, 'r')
    seqs1 = {}
    seqs2 = {}
    seq_id1_dict = {}
    seq_id2_dict = {}
    seq_id = ''
    path = os.getcwd()
    os.mkdir('R_block')
    while True:
        line = inf1.readline()
        if not line: break
        if line[0] == '>':
            seq_id = line[1:].strip()
            seqs1[seq_id] = ''
        else:
            seqs1[seq_id] += line.strip()
    while True:
        line = inf2.readline()
        if not line: break
        if line[0] == '>':
            seq_id = line[1:].strip()
            seqs2[seq_id] = ''
        else:
            seqs2[seq_id] += line.strip()
    for n in synteny_pair:
        pair = synteny_pair[n]
        seq_id1 = pair[1]
        seq_id2 = pair[2]
        outf_query = open(path + '/R_block/' + str(n) + '.query.fa', 'a')
        outf_subject = open(path + '/R_block/' + str(n) + '.subject.fa', 'a')
        outf_query.write('>' + seq_id1 + '\n' + seqs1[seq_id1])
        outf_subject.write('>' + seq_id2 + '\n' + seqs2[seq_id2])
        seq_id1_dict[n] = seq_id1
        seq_id2_dict[n] = seq_id2
        outf_query.close()
        outf_subject.close()
    inf1.close()
    inf2.close()
    nfile = open(path + '/R_block/n.txt', 'a')
    nfile.write(str(num))
    return num, synteny_pair, seq_id1_dict, seq_id2_dict, block_range_sp1, block_range_sp2


def needlecmd(path, n, seq_id1_dict, seq_id2_dict):
    alignlog = open(path + '/R_block/' + str(n) + '.alignlog', 'a')
    needlecommandline = 'needle ' + path + '/R_block/' + str(
        n) + '.query.fa ' + path + '/R_block/' + str(
        n) + '.subject.fa -outfile ' + path + '/R_block/' + str(
        n) + '.needle -gapopen 10.0 -gapextend 0.5 -sprotein1 Y -sprotein2 Y -sid1 ' + \
                        seq_id1_dict[n] + ' -sid2 ' + seq_id2_dict[n]
    alignlog.write('Command:' + needlecommandline + '\n')
    p = subprocess.run(shlex.split(needlecommandline), shell=False, stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT, encoding="utf-8")
    alignlog.write(p.stdout + '\n' + p.stderr)
    alignlog.close()


def needle(inf_anchors, infasta1, infasta2, pnumber, inbed_sp1, inbed_sp2):
    path = os.getcwd()
    n, synteny_pair, seq_id1_dict, seq_id2_dict, block_range_sp1, block_range_sp2 = \
        dict_anchor(inf_anchors, infasta1, infasta2, inbed_sp1, inbed_sp2)
    pool = multiprocessing.Pool(processes=int(pnumber))
    for i in range(int(n)):
        pool.apply_async(needlecmd, (path, i, seq_id1_dict, seq_id2_dict))
    pool.close()
    pool.join()
    return n, synteny_pair, seq_id1_dict, seq_id2_dict, block_range_sp1, block_range_sp2


def reformatneedle(path, n):
    seq_id_dict = {}
    Length_dict = {}
    Identity_dict = {}
    Similarity_dict = {}
    Gaps_dict = {}
    Score_dict = {}
    R_dict = {}
    for i in range(int(n)):
        nf = open(path + '/R_block/' + str(i) + '.needle', 'r')
        seq_id1 = ''
        seq_id2 = ''
        Length = ''
        Identity_counts = ''
        Similarity_counts = ''
        Gaps_counts = ''
        Score = ''
        while True:
            line = nf.readline()
            if not line: break
            if line[0:4] == '# 1:':
                seq_id1 = line[5:].strip('\n')
            elif line[0:4] == '# 2:':
                seq_id2 = line[5:].strip('\n')
            elif line[0:9] == '# Length:':
                Length = line[10:].strip('\n')
            elif line[0:11] == '# Identity:':
                Identity_counts = ((line.split(':')[1]).split('/')[0]).replace(' ', '')
            elif line[0:13] == '# Similarity:':
                Similarity_counts = ((line.split(':')[1]).split('/')[0]).replace(' ', '')
            elif line[0:7] == '# Gaps:':
                Gaps_counts = ((line.split(':')[1]).split('/')[0]).replace(' ', '')
            elif line[0:8] == '# Score:':
                Score = ((line.split(':')[1]).rstrip('\n')).replace(' ', '')
        Identity = int(Identity_counts) / int(Length)
        Similarity = int(Similarity_counts) / int(Length)
        Gaps = int(Gaps_counts) / int(Length)
        R = Similarity * (1 - Gaps)
        seq_id_dict[i] = (seq_id1, seq_id2)
        Length_dict[i] = Length
        Identity_dict[i] = Identity
        Similarity_dict[i] = Similarity
        Gaps_dict[i] = Gaps
        Score_dict[i] = Score
        R_dict[i] = R
        nf.close()
    return seq_id_dict, Length_dict, Identity_dict, Similarity_dict, Gaps_dict, Score_dict, R_dict


def blockR(inf_anchors, sp1blockRbed, sp2blockRbed, infasta1, infasta2, pnumber, inbed_sp1, inbed_sp2, outR):
    path = os.getcwd()
    n, synteny_pair, seq_id1_dict, seq_id2_dict, block_range_sp1, block_range_sp2 = \
        needle(inf_anchors, infasta1, infasta2, pnumber, inbed_sp1, inbed_sp2)
    # numf = open(path + '/R_block/n.txt', 'r')
    # n = numf.readline()
    seq_id_dict, Length_dict, Identity_dict, Similarity_dict, Gaps_dict, Score_dict, R_dict = \
        reformatneedle(path, n)
    blockR_dict = {}
    oR = open(outR, 'w')
    oR.write("Block\tsp1ID\tsp2ID\tLength\tIdentity\tSimilarity\tGaps\tScore\tR\n")
    for i in range(int(n)):
        pair = synteny_pair[i]
        block_name = pair[0]
        seq_id1 = pair[1]
        seq_id2 = pair[2]
        oR.write(block_name + '\t' + seq_id1 + '\t' + seq_id2 + '\t' + Length_dict[i] + '\t' +
                 str(round(Identity_dict[i], 4)) + '\t' + str(round(Similarity_dict[i], 4)) + '\t' +
                 str(round(Gaps_dict[i], 4)) + '\t' + Score_dict[i] + '\t' +
                 str(round(R_dict[i], 4)) + '\n')
        if block_name not in blockR_dict:
            blockR_dict[block_name] = [R_dict[i]]
        elif block_name in blockR_dict:
            blockR_dict[block_name].append(R_dict[i])
    oR.close()
    
    blockR25 = {}
    for block_name in blockR_dict:
        blockR25[block_name] = np.percentile(np.array(blockR_dict[block_name]), 25)
    blockR_dictfilter = {}
    for block_name in blockR25:
        if blockR25[block_name] >= 0.5:
            blockR_dictfilter[block_name] = 0.5
        elif blockR25[block_name] < 0.5:
            blockR_dictfilter[block_name] = blockR25[block_name]
    blockRbed_sp1 = open(sp1blockRbed, 'w')
    blockRbed_sp2 = open(sp2blockRbed, 'w')
    for block_name in block_range_sp1:
        block_region = block_range_sp1[block_name][0] + '\t' + \
                       str(block_range_sp1[block_name][1]) + '\t' + \
                       str(block_range_sp1[block_name][2])
        blockRbed_sp1.write(block_region + '\t' + block_name + '\t' + str(round(blockR_dictfilter[block_name], 4)) + '\n')
    for block_name in block_range_sp2:
        block_region = block_range_sp2[block_name][0] + '\t' + \
                       str(block_range_sp2[block_name][1]) + '\t' + \
                       str(block_range_sp2[block_name][2])
        blockRbed_sp2.write(block_region + '\t' + block_name + '\t' + str(round(blockR_dictfilter[block_name], 4)) + '\n')
    blockRbed_sp1.close()
    blockRbed_sp2.close()


if __name__ == '__main__':
    inf_anchors = sys.argv[1]
    sp1blockRbed = sys.argv[2]
    sp2blockRbed = sys.argv[3]
    infasta1 = sys.argv[4]
    infasta2 = sys.argv[5]
    pnumber = sys.argv[6]
    inbed_sp1 = sys.argv[7]
    inbed_sp2 = sys.argv[8]
    outR = sys.argv[9]
    blockR(inf_anchors, sp1blockRbed, sp2blockRbed, infasta1, infasta2, pnumber, inbed_sp1, inbed_sp2, outR)