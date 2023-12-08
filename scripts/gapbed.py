#!/usr/bin/env python3
# -*- coding:utf-8 -*-
import sys


def dict_bed(infile):
    inf_bed = open(infile, 'r')
    bed_dict1 = {}
    bed_dict2 = {}
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
        if gene_id not in bed_dict1:
            n += 1
            bed_dict1[gene_id] = n
            bed_dict2[n] = gene_id
            gene_chr[gene_id] = chr
            gene_start[gene_id] = start
            gene_end[gene_id] = end
            gene_direction[gene_id] = direction
    inf_bed.close()
    return bed_dict1, bed_dict2, gene_chr, gene_start, gene_end, gene_direction


def gap_dict(infile):
    inf_gap = open(infile, 'r')
    lines = inf_gap.readlines()
    gap_lines = []
    for line in lines:
        if line != '###\n' and line != '#\n':
            gap = line.rstrip('\n')
            gap_lines.append(gap)
        else:
            continue
    return gap_lines


def gapbed(in_gap, inbed_sp1, inbed_sp2, outgapbed_sp1, outgapbed_sp2):
    print('[\033[0;36mINFO\033[0m] Load file `\033[0;35m' + in_gap + '\033[0m`')
    gap_lines = gap_dict(in_gap)
    print('[\033[0;36mINFO\033[0m] Load file `' + inbed_sp1 + '`')
    bed_sp11, bed_sp12, chr_sp1, start_sp1, end_sp1, direction_sp1 = dict_bed(inbed_sp1)
    print('[\033[0;36mINFO\033[0m] Load file `' + inbed_sp2 + '`')
    bed_sp21, bed_sp22, chr_sp2, start_sp2, end_sp2, direction_sp2 = dict_bed(inbed_sp2)
    outf_sp1 = open(outgapbed_sp1, 'w')
    outf_sp2 = open(outgapbed_sp2, 'w')
    gap_sp1 = []
    gap_sp2 = []
    for line in gap_lines:
        gapstart_sp1, gapend_sp1, gapstart_sp2, gapend_sp2 = line.split('\t')
        gap_sp1.append([gapstart_sp1, gapend_sp1])
        gap_sp2.append([gapstart_sp2, gapend_sp2])
    print('[\033[0;36mINFO\033[0m] Extracting ...')
    for n in range(len(gap_sp1)):
        gap_1n = gap_sp1[n]
        gap_2n = gap_sp2[n]
        gap1_gene_number = abs(bed_sp11[gap_1n[1]] - bed_sp11[gap_1n[0]])
        if int(bed_sp11[gap_1n[0]]) < int(bed_sp11[gap_1n[1]]):
            for m in range(gap1_gene_number + 1):
                sp1m = bed_sp11[gap_1n[0]] + m
                if sp1m < bed_sp11[gap_1n[1]]:
                    outf_sp1.write('1' + '\t' + chr_sp1[bed_sp12[sp1m]] + '\t' + start_sp1[bed_sp12[sp1m]] + '\t' +
                                   end_sp1[bed_sp12[sp1m]] + '\t' + bed_sp12[sp1m] + ';')
                elif sp1m == bed_sp11[gap_1n[1]]:
                    outf_sp1.write('1' + '\t' + chr_sp1[bed_sp12[sp1m]] + '\t' + start_sp1[bed_sp12[sp1m]] + '\t' +
                                   end_sp1[bed_sp12[sp1m]] + '\t' + bed_sp12[sp1m] + '\n')
            if int(start_sp1[gap_1n[1]]) >= int(end_sp1[gap_1n[0]]):
                outf_sp1.write('2' + '\t' + chr_sp1[gap_1n[0]] + '\t' + end_sp1[gap_1n[0]] + '\t' +
                               start_sp1[gap_1n[1]] + '\t' + gap_1n[0] + '_' + gap_1n[1] + '\n')
            elif int(end_sp1[gap_1n[0]]) > int(start_sp1[gap_1n[1]]):
                start = []
                end = []
                for m in range(gap1_gene_number):
                    sp1m = bed_sp11[gap_1n[0]] + m
                    start.append(int(start_sp1[bed_sp12[sp1m]]))
                    end.append(int(end_sp1[bed_sp12[sp1m]]))
                outf_sp1.write('2' + '\t' + chr_sp1[gap_1n[0]] + '\t' + str(min(start)) + '\t' + str(max(end)) + '\t' +
                               gap_1n[0] + '_' + gap_1n[1] + '\n')
        elif bed_sp11[gap_1n[0]] > bed_sp11[gap_1n[1]]:
            for m in range(gap1_gene_number + 1):
                sp1m = bed_sp11[gap_1n[1]] + m
                if sp1m < bed_sp11[gap_1n[0]]:
                    outf_sp1.write('1' + '\t' + chr_sp1[bed_sp12[sp1m]] + '\t' + start_sp1[bed_sp12[sp1m]] + '\t' +
                                   end_sp1[bed_sp12[sp1m]] + '\t' + bed_sp12[sp1m] + ';')
                elif sp1m == bed_sp11[gap_1n[0]]:
                    outf_sp1.write('1' + '\t' + chr_sp1[bed_sp12[sp1m]] + '\t' + start_sp1[bed_sp12[sp1m]] + '\t' +
                                   end_sp1[bed_sp12[sp1m]] + '\t' + bed_sp12[sp1m] + '\n')
            if int(start_sp1[gap_1n[0]]) >= int(end_sp1[gap_1n[1]]):
                outf_sp1.write('2' + '\t' + chr_sp1[gap_1n[0]] + '\t' + end_sp1[gap_1n[1]] + '\t' +
                               start_sp1[gap_1n[0]] + '\t' + gap_1n[1] + '_' + gap_1n[0] + '\n')
            elif int(end_sp1[gap_1n[1]]) > int(start_sp1[gap_1n[0]]):
                start = []
                end = []
                for m in range(gap1_gene_number):
                    sp1m = bed_sp11[gap_1n[1]] + m
                    start.append(int(start_sp1[bed_sp12[sp1m]]))
                    end.append(int(end_sp1[bed_sp12[sp1m]]))
                outf_sp1.write('2' + '\t' + chr_sp1[gap_1n[0]] + '\t' + str(min(start)) + '\t' + str(max(end)) + '\t' +
                               gap_1n[1] + '_' + gap_1n[0] + '\n')
        gap2_gene_number = abs(bed_sp21[gap_2n[1]] - bed_sp21[gap_2n[0]])
        if bed_sp21[gap_2n[0]] < bed_sp21[gap_2n[1]]:
            if int(start_sp2[gap_2n[1]]) >= int(end_sp2[gap_2n[0]]):
                outf_sp2.write('1' + '\t' + chr_sp2[gap_2n[0]] + '\t' + end_sp2[gap_2n[0]] + '\t' +
                               start_sp2[gap_2n[1]] + '\t' + gap_2n[0] + '_' + gap_2n[1] + '\n')
            elif int(end_sp2[gap_2n[0]]) > int(start_sp2[gap_2n[1]]):
                start = []
                end = []
                for m in range(gap2_gene_number):
                    sp2m = bed_sp21[gap_2n[0]] + m
                    start.append(int(start_sp2[bed_sp22[sp2m]]))
                    end.append(int(end_sp2[bed_sp22[sp2m]]))
                outf_sp2.write('1' + '\t' + chr_sp2[gap_2n[0]] + '\t' + str(min(start)) + '\t' + str(max(end)) + '\t' +
                               gap_2n[0] + '_' + gap_2n[1] + '\n')
            for m in range(gap2_gene_number + 1):
                sp2m = bed_sp21[gap_2n[0]] + m
                if sp2m < bed_sp21[gap_2n[1]]:
                    outf_sp2.write('2' + '\t' + chr_sp2[bed_sp22[sp2m]] + '\t' + start_sp2[bed_sp22[sp2m]]
                                   + '\t' + end_sp2[bed_sp22[sp2m]] + '\t' + bed_sp22[sp2m] + ';')
                elif sp2m == bed_sp21[gap_2n[1]]:
                    outf_sp2.write('2' + '\t' + chr_sp2[bed_sp22[sp2m]] + '\t' + start_sp2[bed_sp22[sp2m]]
                                   + '\t' + end_sp2[bed_sp22[sp2m]] + '\t' + bed_sp22[sp2m] + '\n')
        elif bed_sp21[gap_2n[0]] > bed_sp21[gap_2n[1]]:
            if int(start_sp2[gap_2n[0]]) >= int(end_sp2[gap_2n[1]]):
                outf_sp2.write('1' + '\t' + chr_sp2[gap_2n[0]] + '\t' + end_sp2[gap_2n[1]] + '\t' +
                               start_sp2[gap_2n[0]] + '\t' + gap_2n[1] + '_' + gap_2n[0] + '\n')
            elif int(end_sp2[gap_2n[1]]) > int(start_sp2[gap_2n[0]]):
                start = []
                end = []
                for m in range(gap2_gene_number):
                    sp2m = bed_sp21[gap_2n[0]] + m
                    start.append(int(start_sp2[bed_sp22[sp2m]]))
                    end.append(int(end_sp2[bed_sp22[sp2m]]))
                outf_sp2.write('1' + '\t' + chr_sp2[gap_2n[0]] + '\t' + str(min(start)) + '\t' + str(max(end)) + '\t' +
                               gap_2n[1] + '_' + gap_2n[1] + '\n')
            for m in range(gap2_gene_number + 1):
                sp2m = bed_sp21[gap_2n[1]] + m
                if sp2m < bed_sp21[gap_2n[0]]:
                    outf_sp2.write('2' + '\t' + chr_sp2[bed_sp22[sp2m]] + '\t' + start_sp2[bed_sp22[sp2m]]
                                   + '\t' + end_sp2[bed_sp22[sp2m]] + '\t' + bed_sp22[sp2m] + ';')
                elif sp2m == bed_sp21[gap_2n[0]]:
                    outf_sp2.write('2' + '\t' + chr_sp2[bed_sp22[sp2m]] + '\t' + start_sp2[bed_sp22[sp2m]]
                                   + '\t' + end_sp2[bed_sp22[sp2m]] + '\t' + bed_sp22[sp2m] + '\n')
        elif bed_sp21[gap_2n[0]] == bed_sp21[gap_2n[1]]:
            outf_sp2.write('1' + '\t' + chr_sp2[gap_2n[0]] + '\t' + start_sp2[gap_2n[0]] + '\t' +
                           end_sp2[gap_2n[0]] + '\t' + gap_2n[0] + '_' + gap_2n[0] + '\n')
            outf_sp2.write('1' + '\t' + chr_sp2[gap_2n[0]] + '\t' + start_sp2[gap_2n[0]] + '\t' +
                           end_sp2[gap_2n[0]] + '\t' + gap_2n[0] + '\n')
    outf_sp1.close()
    outf_sp2.close()


if __name__ == '__main__':
    in_gap = sys.argv[1]
    inbed_sp1 = sys.argv[2]
    inbed_sp2 = sys.argv[3]
    outgapbed_sp1 = sys.argv[4]
    outgapbed_sp2 = sys.argv[5]
    gapbed(in_gap, inbed_sp1, inbed_sp2, outgapbed_sp1, outgapbed_sp2)
