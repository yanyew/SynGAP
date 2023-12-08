#!/usr/bin/env python3
# -*- coding:utf-8 -*-
import sys
from collections import Counter


def dict_bed(infile):
    inf_bed = open(infile, 'r')
    bed_dict = {}
    gene_chr = {}
    gene_start = {}
    gene_end = {}
    gene_direction = {}
    n = -1
    while True:
        line = inf_bed.readline()
        if not line: break
        chr = line.split('\t')[0]
        start = line.split('\t')[1]
        end = line.split('\t')[2]
        gene_id = line.split('\t')[3]
        direction = line.split('\t')[5]
        if gene_id not in bed_dict:
            n += 1
            bed_dict[gene_id] = n
            gene_chr[gene_id] = chr
            gene_start[gene_id] = start
            gene_end[gene_id] = end
            gene_direction[gene_id] = direction
    inf_bed.close()
    return bed_dict, gene_chr, gene_start, gene_end, gene_direction


def list_blocks(inf_anchors):
    in_anchors = open(inf_anchors, 'r')
    anchors = in_anchors.read()
    anchors_list = anchors.split('###\n')
    blocks_list = []
    for block in anchors_list:
        if block != '':
            blocks_list.append(block.rstrip('\n'))
        else:
            continue
    in_anchors.close()
    return blocks_list


def dict_block(block, inbed_sp1, inbed_sp2):
    bed_sp1, chr_sp1, start_sp1, end_sp1, direction_sp1 = dict_bed(inbed_sp1)
    bed_sp2, chr_sp2, start_sp2, end_sp2, direction_sp2 = dict_bed(inbed_sp2)
    keys_bed_sp1 = list(bed_sp1.keys())
    keys_bed_sp2 = list(bed_sp2.keys())
    sp1_geneid_list = []
    block_dict_list = []
    block_sp1_geneid = []
    block_sp1_geneid_list = []
    block_line_sp1 = {}
    block_line_sp1_list = []
    block_line_sp2 = {}
    block_line_sp2_list = []
    line_num = 0
    repeat = 1
    block_derection = 1
    lines = block.split('\n')
    while '' in lines:
        lines.remove('')
    for line in lines:
        species_sp1_gene_id = line.split('\t')[0]
        block_line_sp1[line_num] = species_sp1_gene_id
        sp1_geneid_list.append(species_sp1_gene_id)
        species_sp2_gene_id = line.split('\t')[1]
        block_line_sp2[line_num] = species_sp2_gene_id
        line_num += 1
    keys_block_line_sp1 = list(block_line_sp1.keys())
    values_block_line_sp1 = list(block_line_sp1.values())
    values_block_line_sp2 = list(block_line_sp2.values())
    repeat_times = Counter(sp1_geneid_list)
    repeat_times_dict = dict(repeat_times)
    max_repeat_key = repeat_times.most_common(1)
    max_repeat = max_repeat_key[0][1]
    for m in range(max_repeat):
        x = locals()['block_line_sp1_' + str(m)] = {}
        block_line_sp1_list.append(x.copy())
        y = locals()['block_line_sp2_' + str(m)] = {}
        block_line_sp2_list.append(y.copy())
        z = locals()['block_dict_' + str(m)] = {}
        block_dict_list.append(z.copy())
        u = locals()['block_sp1_geneid_' + str(m)] = []
        block_sp1_geneid_list.append(u.copy())
    for key in repeat_times_dict:
        for n in range(max_repeat):
            if (n + 1) <= repeat_times_dict[key]:
                block_sp1_geneid_list[n].append(key)
            else:
                break
    for n in range(len(block_line_sp1)):
        sp1_geneid = block_line_sp1[n]
        sp2_geneid = block_line_sp2[n]
        if repeat_times_dict[block_line_sp1[n]] == 1:
            repeat = 1
            position_0 = block_sp1_geneid_list[0].index(sp1_geneid)
            if (n - 1) < 0:
                block_line_sp1_list[0][position_0] = sp1_geneid
                block_line_sp2_list[0][position_0] = sp2_geneid
                block_dict_list[0][sp1_geneid] = sp2_geneid
            elif (n - 1) >= 0:
                sp1_geneid_p = block_sp1_geneid_list[0][position_0 - 1]
                if repeat_times_dict[sp1_geneid_p] == 1:
                    block_line_sp1_list[0][position_0] = sp1_geneid
                    block_line_sp2_list[0][position_0] = sp2_geneid
                    block_dict_list[0][sp1_geneid] = sp2_geneid
                elif repeat_times_dict[sp1_geneid_p] > 1:
                    sp2_geneid_p = block_dict_list[0][sp1_geneid_p]
                    sp2_geneid_p1 = block_dict_list[repeat][sp1_geneid_p]
                    bed_position_sp2 = bed_sp2[sp2_geneid]
                    bed_position_sp2_p = bed_sp2[sp2_geneid_p]
                    bed_position_sp2_p1 = bed_sp2[sp2_geneid_p1]
                    if abs(bed_position_sp2 - bed_position_sp2_p) < abs(
                            bed_position_sp2 - bed_position_sp2_p1):
                        block_line_sp1_list[0][position_0] = sp1_geneid
                        block_line_sp2_list[0][position_0] = sp2_geneid
                        block_dict_list[0][sp1_geneid] = sp2_geneid
                    elif abs(bed_position_sp2 - bed_position_sp2_p) >= abs(
                            bed_position_sp2 - bed_position_sp2_p1):
                        block_sp1_geneid_list[0].remove(sp1_geneid)
                        block_sp1_geneid_list[repeat].append(sp1_geneid)
                        position_repeat = block_sp1_geneid_list[repeat].index(sp1_geneid)
                        block_line_sp1_list[repeat][position_repeat] = sp1_geneid
                        block_line_sp2_list[repeat][position_repeat] = sp2_geneid
                        block_dict_list[repeat][sp1_geneid] = sp2_geneid
        elif repeat_times_dict[block_line_sp1[n]] > 1:
            position_0 = block_sp1_geneid_list[0].index(sp1_geneid)
            if (n - 1) < 0:
                repeat = 1
                block_line_sp1_list[0][position_0] = sp1_geneid
                block_line_sp2_list[0][position_0] = sp2_geneid
                block_dict_list[0][sp1_geneid] = sp2_geneid
                repeat += 1
            elif (n - 1) >= 0:
                if block_line_sp1[n] != block_line_sp1[n - 1]:
                    repeat = 1
                    position_repeat_1 = block_sp1_geneid_list[repeat].index(sp1_geneid)
                    block_position = values_block_line_sp1.index(sp1_geneid)
                    sp1_geneid_p = block_sp1_geneid_list[0][position_0 - 1]
                    sp2_geneid_p = block_dict_list[0][sp1_geneid_p]
                    sp2_geneid_1 = block_line_sp2[block_position + 1]
                    bed_position_sp2 = bed_sp2[sp2_geneid]
                    bed_position_sp2_p = bed_sp2[sp2_geneid_p]
                    bed_position_sp2_1 = bed_sp2[sp2_geneid_1]
                    if abs(bed_position_sp2 - bed_position_sp2_p) < abs(
                            bed_position_sp2_1 - bed_position_sp2_p):
                        block_line_sp1_list[0][position_0] = sp1_geneid
                        block_line_sp2_list[0][position_0] = sp2_geneid
                        block_dict_list[0][sp1_geneid] = sp2_geneid
                        repeat += 1
                        block_derection = 1
                    elif abs(bed_position_sp2 - bed_position_sp2_p) >= abs(
                            bed_position_sp2_1 - bed_position_sp2_p):
                        block_line_sp1_list[0][position_0] = sp1_geneid
                        block_line_sp2_list[0][position_0] = sp2_geneid
                        block_dict_list[0][sp1_geneid] = sp2_geneid
                        repeat += 1
                        block_derection = -1
                elif block_line_sp1[n] == block_line_sp1[n - 1]:
                    position_repeat_1 = block_sp1_geneid_list[repeat - 1].index(sp1_geneid)
                    block_position = values_block_line_sp1.index(sp1_geneid) + repeat - 1
                    if position_0 > 0:
                        if block_derection == 1:
                            block_line_sp1_list[repeat - 1][position_repeat_1] = sp1_geneid
                            block_line_sp2_list[repeat - 1][position_repeat_1] = sp2_geneid
                            block_dict_list[repeat - 1][sp1_geneid] = sp2_geneid
                            repeat += 1
                        elif block_derection == -1:
                            sp1_geneid_p = block_sp1_geneid_list[0][position_0 - 1]
                            sp2_geneid_p = block_dict_list[0][sp1_geneid_p]
                            sp2_geneid_0 = block_dict_list[0][sp1_geneid]
                            bed_position_sp2 = bed_sp2[sp2_geneid]
                            bed_position_sp2_p = bed_sp2[sp2_geneid_p]
                            bed_position_sp2_0 = bed_sp2[sp2_geneid_0]
                            if abs(bed_position_sp2 - bed_position_sp2_p) < abs(
                                    bed_position_sp2_0 - bed_position_sp2_p):
                                block_line_sp1_list[0][position_0] = sp1_geneid
                                block_line_sp2_list[0][position_0] = sp2_geneid
                                block_dict_list[0][sp1_geneid] = sp2_geneid
                                block_line_sp1_list[repeat - 1][position_repeat_1] = sp1_geneid
                                block_line_sp2_list[repeat - 1][position_repeat_1] = sp2_geneid_0
                                block_dict_list[repeat - 1][sp1_geneid] = sp2_geneid_0
                                repeat += 1
                            elif abs(bed_position_sp2 - bed_position_sp2_p) >= abs(
                                    bed_position_sp2_0 - bed_position_sp2_p):
                                block_line_sp1_list[repeat - 1][position_repeat_1] = sp1_geneid
                                block_line_sp2_list[repeat - 1][position_repeat_1] = sp2_geneid
                                block_dict_list[repeat - 1][sp1_geneid] = sp2_geneid
                                repeat += 1
                    elif position_0 == 0:
                        block_line_sp1_list[repeat - 1][position_repeat_1] = sp1_geneid
                        block_line_sp2_list[repeat - 1][position_repeat_1] = sp2_geneid
                        block_dict_list[repeat - 1][sp1_geneid] = sp2_geneid
                        repeat += 1
    return block_line_sp1_list, block_line_sp2_list, block_dict_list, max_repeat


def gap_find(inbed_sp1, inbed_sp2, inf_anchors, outfile):
    outf = open(outfile, 'w')
    print('[\033[0;36mINFO\033[0m] Load file `\033[0;35m' + inf_anchors + '\033[0m`')
    blocks_list = list_blocks(inf_anchors)
    print('[\033[0;36mINFO\033[0m] Load file `' + inbed_sp1 + '`')
    bed_sp1, chr_sp1, start_sp1, end_sp1, direction_sp1 = dict_bed(inbed_sp1)
    print('[\033[0;36mINFO\033[0m] Load file `' + inbed_sp2 + '`')
    bed_sp2, chr_sp2, start_sp2, end_sp2, direction_sp2 = dict_bed(inbed_sp2)
    for block in blocks_list:
        block_line_sp1_list, block_line_sp2_list, block_dict_list, max_repeat = \
            dict_block(block, inbed_sp1, inbed_sp2)
        num = 1
        for block_line_sp1_n, block_line_sp2_n, block_dict_n in zip(
                block_line_sp1_list, block_line_sp2_list, block_dict_list):
            keys_block_line_sp1_n = list(block_line_sp1_n.keys())
            for key in block_line_sp1_n:
                sp1_geneid = block_line_sp1_n[key]
                sp2_geneid = block_line_sp2_n[key]
                position = keys_block_line_sp1_n.index(key)
                if position in range(len(block_line_sp1_n) - 1):
                    key1 = keys_block_line_sp1_n[position + 1]
                    sp1_geneid1 = block_line_sp1_n[key1]
                    sp2_geneid1 = block_line_sp2_n[key1]
                    if (abs(bed_sp1[sp1_geneid1] - bed_sp1[sp1_geneid]) != 1 and
                        abs(bed_sp1[sp1_geneid1] - bed_sp1[sp1_geneid]) <= 20 and
                        abs(bed_sp2[sp2_geneid1] - bed_sp2[sp2_geneid]) <= 20) or \
                            (abs(bed_sp2[sp2_geneid1] - bed_sp2[sp2_geneid]) != 1 and
                             abs(bed_sp1[sp1_geneid1] - bed_sp1[sp1_geneid]) <= 20 and
                             abs(bed_sp2[sp2_geneid1] - bed_sp2[sp2_geneid]) <= 20):
                        outf.write(sp1_geneid + '\t' + sp1_geneid1 + '\t' +
                                   block_dict_n[sp1_geneid] + '\t' + block_dict_n[sp1_geneid1] + '\n')
            if num < max_repeat:
                outf.write('#' + '\n')
                num += 1
            elif num == max_repeat:
                outf.write('###' + '\n')
    outf.close()


if __name__ == '__main__':
    inbed_sp1 = sys.argv[1]
    inbed_sp2 = sys.argv[2]
    inf_anchors = sys.argv[3]
    outfile = sys.argv[4]
    gap_find(inbed_sp1, inbed_sp2, inf_anchors, outfile)
