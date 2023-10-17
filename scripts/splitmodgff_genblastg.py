#!/usr/bin/env python3
# -*- coding:utf-8 -*-
import sys


def dict_gapbed(infile):
    inf_gapbed = open(infile, 'r')
    gap_chr = {}
    gap_start = {}
    gap_end = {}
    while True:
        line = inf_gapbed.readline()
        if not line: break
        chr = line.split('\t')[1]
        start = line.split('\t')[2]
        end = line.split('\t')[3]
        gap_id = line.split('\t')[4].rstrip('\n')
        if '_' in gap_id:
            gap_chr[gap_id] = chr
            gap_start[gap_id] = start
            gap_end[gap_id] = end
    inf_gapbed.close()
    return gap_chr, gap_start, gap_end


def split(sp1, sp2, ingapbed_sp1, ingapbed_sp2, inmodgff, outgff_sp1, outgff_sp2):
    outfgff_sp1 = open(outgff_sp1, 'w')
    outfgff_sp2 = open(outgff_sp2, 'w')
    print('[\033[0;36mINFO\033[0m] Load file `\033[0;35m' + inmodgff + '\033[0m`')
    inf_gff = open(inmodgff, 'r')
    lines = inf_gff.readlines()
    print('[\033[0;36mINFO\033[0m] Load file `' + ingapbed_sp1 + '`')
    gap_chr_1, gap_start_1, gap_end_1 = dict_gapbed(ingapbed_sp1)
    print('[\033[0;36mINFO\033[0m] Load file `' + ingapbed_sp2 + '`')
    gap_chr_2, gap_start_2, gap_end_2 = dict_gapbed(ingapbed_sp2)
    name_list = []
    transcripts = {}
    rename = {}
    print('[\033[0;36mINFO\033[0m] Spliting ...')
    for line in lines:
        if line[0] != '#':
            type = line.split('\t')[2]
            attributes = (line.split('\t')[8]).rstrip('\n')
            if type == 'transcript':
                ID = (attributes.split(';')[0]).split('=')[1]
                Name = (attributes.split(';')[1]).split('=')[1]
                rename[ID] = Name
    for line in lines:
        if line[0] != '#':
            type = line.split('\t')[2]
            attributes = (line.split('\t')[8]).rstrip('\n')
            if type == 'transcript':
                ID = (attributes.split(';')[0]).split('=')[1]
                Name = rename[ID]
                name_list.append(Name)
                if name_list.count(Name) == 1:
                    Name_alter = Name + '_1'
                    line_alter = line.replace(ID, Name_alter)
                    line_alter_1 = line_alter.split(';Name=')[0] + ';Name=' + Name_alter + ';Target=' + Name + '\n'
                    transcripts[Name_alter] = ''
                    transcripts[Name_alter] += line_alter_1
                elif name_list.count(Name) > 1:
                    n = name_list.count(Name)
                    Name_alter = Name + '_' + str(n)
                    line_alter = line.replace(ID, Name_alter)
                    line_alter_1 = line_alter.split(';Name=')[0] + ';Name=' + Name_alter + ';Target=' + Name + '\n'
                    transcripts[Name_alter] = ''
                    transcripts[Name_alter] += line_alter_1
            elif type == 'coding_exon':
                Parent = (attributes.split(';')[1]).split('=')[1]
                Name = rename[Parent]
                # line_rename = line.replace(Parent, Name)
                if name_list.count(Name) == 1:
                    Name_alter = Name + '_1'
                    line_alter = line.replace(Parent, Name_alter)
                    line_alter_1 = line_alter.split(';Parent=')[0] + ';Parent=' + Name_alter + ';Target=' + Name + '\n'
                    transcripts[Name_alter] += line_alter_1
                elif name_list.count(Name) > 1:
                    n = name_list.count(Name)
                    Name_alter = Name + '_' + str(n)
                    line_alter = line.replace(Parent, Name_alter)
                    line_alter_1 = line_alter.split(';Parent=')[0] + ';Parent=' + Name_alter + ';Target=' + Name + '\n'
                    transcripts[Name_alter] += line_alter_1
    for name in transcripts:
        for line in (transcripts[name]).split('\n'):
            if line != '':
                gapid = line.split('\t')[0]
                source = line.split('\t')[1]
                type = line.split('\t')[2]
                start = line.split('\t')[3]
                end = line.split('\t')[4]
                score = line.split('\t')[5]
                strand = line.split('\t')[6]
                phase = line.split('\t')[7]
                attributes = (line.split('\t')[8]).rstrip('\n')
                if type == 'transcript':
                    type = 'mRNA'
                    if gapid in gap_chr_1:
                        attributes = attributes.replace('ID=', 'ID=' + sp1 + '_')
                        attributes = attributes.replace('Name=', 'Name=' + sp1 + '_')
                        outfgff_sp1.write(gap_chr_1[gapid] + '\t' + source + '\t' + type + '\t' + str(
                            int(start) + int(gap_start_1[gapid]))
                                          + '\t' + str(
                            int(end) + int(gap_start_1[gapid])) + '\t' + score + '\t' + strand + '\t' + phase
                                          + '\t' + attributes + '\n')
                    elif gapid in gap_chr_2:
                        attributes = attributes.replace('ID=', 'ID=' + sp2 + '_')
                        attributes = attributes.replace('Name=', 'Name=' + sp2 + '_')
                        outfgff_sp2.write(gap_chr_2[gapid] + '\t' + source + '\t' + type + '\t' + str(
                            int(start) + int(gap_start_2[gapid]))
                                             + '\t' + str(
                            int(end) + int(gap_start_2[gapid])) + '\t' + score + '\t' + strand + '\t' + phase
                                             + '\t' + attributes + '\n')
                elif type == 'coding_exon':
                    type = 'CDS'
                    if gapid in gap_chr_1:
                        attributes = attributes.replace('ID=', 'ID=' + sp1 + '_')
                        attributes = attributes.replace('Parent=', 'Parent=' + sp1 + '_')
                        outfgff_sp1.write(gap_chr_1[gapid] + '\t' + source + '\t' + type + '\t' + str(
                            int(start) + int(gap_start_1[gapid]))
                                          + '\t' + str(
                            int(end) + int(gap_start_1[gapid])) + '\t' + score + '\t' + strand + '\t' + phase
                                          + '\t' + attributes + '\n')
                    elif gapid in gap_chr_2:
                        attributes = attributes.replace('ID=', 'ID=' + sp2 + '_')
                        attributes = attributes.replace('Parent=', 'Parent=' + sp2 + '_')
                        outfgff_sp2.write(gap_chr_2[gapid] + '\t' + source + '\t' + type + '\t' + str(
                            int(start) + int(gap_start_2[gapid]))
                                             + '\t' + str(
                            int(end) + int(gap_start_2[gapid])) + '\t' + score + '\t' + strand + '\t' + phase
                                             + '\t' + attributes + '\n')
    outfgff_sp1.close()
    outfgff_sp2.close()


if __name__ == '__main__':
    sp1 = sys.argv[1]
    sp2 = sys.argv[2]
    ingapbed_sp1 = sys.argv[3]
    ingapbed_sp2 = sys.argv[4]
    inmodgff = sys.argv[5]
    outgff_sp1 = sys.argv[6]
    outgff_sp2 = sys.argv[7]
    split(sp1, sp2, ingapbed_sp1, ingapbed_sp2, inmodgff, outgff_sp1, outgff_sp2)