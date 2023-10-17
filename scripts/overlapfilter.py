#!/usr/bin/env python3
# -*- coding:utf-8 -*-
import sys


def dict_genomegff(in_genomegff, annotype, annokey):
    inf_genomegff = open(in_genomegff, 'r')
    gene_chrid = {}
    gene_start = {}
    gene_end = {}
    gene_cdslength = {}
    while True:
        line = inf_genomegff.readline()
        if not line: break
        if line[0] != '#' and len(line.split('\t')) == 9:
            chrid = line.split('\t')[0]
            type = line.split('\t')[2]
            start = line.split('\t')[3]
            end = line.split('\t')[4]
            attributes = (line.split('\t')[8]).rstrip('\n')
            attributes_list = attributes.split(';')
            attributes_dict = {}
            for a in attributes_list:
                if a != '':
                    attributes_dict[a.split('=')[0]] = a.split('=')[1]
            if type == str(annotype):
                name = attributes_dict[annokey]
                gene_chrid[name] = chrid
                gene_start[name] = start
                gene_end[name] = end
                if name not in gene_cdslength:
                    gene_cdslength[name] = 0
            elif type == 'CDS':
                Parent = attributes_dict['Parent']
                if Parent not in gene_cdslength:
                    gene_cdslength[Parent] = 0
                    gene_cdslength[Parent] += int(end) - int(start)
                elif Parent in gene_cdslength:
                    gene_cdslength[Parent] += int(end) - int(start)
    inf_genomegff.close()
    return gene_chrid, gene_start, gene_end, gene_cdslength


def dict_modifiedgff(infile):
    inf_modifiedgff = open(infile, 'r')
    gene_chrid = {}
    gene_source = {}
    gene_type = {}
    gene_start = {}
    gene_end = {}
    gene_score = {}
    gene_strand = {}
    gene_phase = {}
    gene_name = {}
    transcripts = []
    mRNAs = {}
    while True:
        line = inf_modifiedgff.readline()
        if not line: break
        chrid = line.split('\t')[0]
        sourse = line.split('\t')[1]
        type = line.split('\t')[2]
        start = line.split('\t')[3]
        end = line.split('\t')[4]
        score = line.split('\t')[5]
        strand = line.split('\t')[6]
        phase = line.split('\t')[7]
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
            gene_chrid[attributes] = chrid
            gene_source[attributes] = sourse
            gene_type[attributes] = type
            gene_start[attributes] = start
            gene_end[attributes] = end
            gene_score[attributes] = score
            gene_strand[attributes] = strand
            gene_phase[attributes] = phase
            gene_name[attributes] = name
        elif type == 'CDS':
            name = attributes_dict['Parent']
            if name in mRNAs:
                mRNAs[name] += line
            elif name not in mRNAs:
                mRNAs[name] = ''
                mRNAs[name] += line
            gene_chrid[attributes] = chrid
            gene_source[attributes] = sourse
            gene_type[attributes] = type
            gene_start[attributes] = start
            gene_end[attributes] = end
            gene_score[attributes] = score
            gene_strand[attributes] = strand
            gene_phase[attributes] = phase
            gene_name[attributes] = name
    inf_modifiedgff.close()
    return gene_chrid, gene_source, gene_type, gene_start, gene_end, gene_score, gene_strand, \
           gene_phase, gene_name, transcripts, mRNAs


def overlapfilter(genomegff, annotype, annokey, modified_gff, filtered_gff):
    print('[\033[0;36mINFO\033[0m] Load file `\033[0;35m' + genomegff + '\033[0m`')
    genome_chrid, genome_start, genome_end, genome_cdslength = dict_genomegff(genomegff, annotype, annokey)
    print('[\033[0;36mINFO\033[0m] Load file `\033[0;35m' + modified_gff + '\033[0m`')
    gene_chrid, gene_source, gene_type, gene_start, gene_end, gene_score, gene_strand, gene_phase, gene_name,\
    transcripts, mRNAs = dict_modifiedgff(modified_gff)
    outfgff = open(filtered_gff, 'w')
    overlap_mRNAs = []
    remain_overlap_mRNAs = []
    overlap_mRNA_length = {}
    print('[\033[0;36mINFO\033[0m] Filtering ...')
    for attributes in gene_chrid:
        if 'mRNA' in gene_type[attributes]:
            for key in genome_chrid:
                if gene_chrid[attributes] == genome_chrid[key]:
                    if int(gene_start[attributes]) <= int(genome_start[key]) and \
                            int(genome_start[key]) <= int(gene_end[attributes]):
                        if gene_name[attributes] in transcripts:
                            transcripts.remove(gene_name[attributes])
                            overlap_mRNAs.append(attributes)
                    elif int(gene_start[attributes]) <= int(genome_end[key]) and \
                            int(genome_end[key]) <= int(gene_end[attributes]):
                        if gene_name[attributes] in transcripts:
                            transcripts.remove(gene_name[attributes])
                            overlap_mRNAs.append(attributes)
                    elif int(gene_start[attributes]) >= int(genome_start[key]) and \
                            int(genome_end[key]) >= int(gene_end[attributes]):
                        if gene_name[attributes] in transcripts:
                            transcripts.remove(gene_name[attributes])
                            overlap_mRNAs.append(attributes)
        else:
            continue
    for attributes in overlap_mRNAs:
        name = (attributes.split(';')[0]).split('=')[1]
        overlap_mRNA_length[name] = 0
        overlap_block = []
        counts = 0
        for line in (mRNAs[name]).split('\n'):
            if line != '':
                type = line.split('\t')[2]
                start = line.split('\t')[3]
                end = line.split('\t')[4]
                if type == 'CDS':
                    length = int(end) - int(start)
                    overlap_mRNA_length[name] += length
        for key in genome_chrid:
            if gene_chrid[attributes] == genome_chrid[key]:
                if int(gene_start[attributes]) <= int(genome_start[key]) and \
                        int(genome_start[key]) <= int(gene_end[attributes]):
                    overlap_block.append(key)
                elif int(gene_start[attributes]) <= int(genome_end[key]) and \
                        int(genome_end[key]) <= int(gene_end[attributes]):
                    overlap_block.append(key)
                elif int(gene_start[attributes]) >= int(genome_start[key]) and \
                        int(genome_end[key]) >= int(gene_end[attributes]):
                    overlap_block.append(key)
        for key in overlap_block:
            if overlap_mRNA_length[name] > (1.1 * genome_cdslength[key]):
                counts += 1
        if counts == len(overlap_block) and gene_name[attributes] not in remain_overlap_mRNAs:
            remain_overlap_mRNAs.append(gene_name[attributes])
    for attributes in gene_chrid:
        attributes_list = attributes.split(';')
        attributes_dict = {}
        for a in attributes_list:
            if a != '':
                attributes_dict[a.split('=')[0]] = a.split('=')[1]
        for transcript in transcripts:
            if transcript not in genome_chrid.keys():
                if gene_type[attributes] == 'mRNA':
                    if transcript == attributes_dict['ID']:
                        outfgff.write(gene_chrid[attributes] + '\t' + gene_source[attributes] + '\t' + gene_type[attributes]
                                  + '\t' + gene_start[attributes] + '\t' + gene_end[attributes] + '\t'
                                  + gene_score[attributes] + '\t' + gene_strand[attributes] + '\t' + gene_phase[attributes]
                                  + '\t' + attributes + ';polish_type=miss_annotated\n')
                elif gene_type[attributes] == 'CDS':
                    if transcript == attributes_dict['Parent']:
                        outfgff.write(gene_chrid[attributes] + '\t' + gene_source[attributes] + '\t' + gene_type[attributes]
                                  + '\t' + gene_start[attributes] + '\t' + gene_end[attributes] + '\t'
                                  + gene_score[attributes] + '\t' + gene_strand[attributes] + '\t' + gene_phase[attributes]
                                  + '\t' + attributes + '\n')
        for mRNA in remain_overlap_mRNAs:
            if mRNA not in genome_chrid.keys():
                if gene_type[attributes] == 'mRNA':
                    if mRNA == attributes_dict['ID']:
                        outfgff.write(gene_chrid[attributes] + '\t' + gene_source[attributes] + '\t' + gene_type[attributes]
                                  + '\t' + gene_start[attributes] + '\t' + gene_end[attributes] + '\t'
                                  + gene_score[attributes] + '\t' + gene_strand[attributes] + '\t' + gene_phase[attributes]
                                  + '\t' + attributes + ';polish_type=mis_annotated\n')
                elif gene_type[attributes] == 'CDS':
                    if mRNA == attributes_dict['Parent']:
                        outfgff.write(gene_chrid[attributes] + '\t' + gene_source[attributes] + '\t' + gene_type[attributes]
                                  + '\t' + gene_start[attributes] + '\t' + gene_end[attributes] + '\t'
                                  + gene_score[attributes] + '\t' + gene_strand[attributes] + '\t' + gene_phase[attributes]
                                  + '\t' + attributes + '\n')
    outfgff.close()


if __name__ == '__main__':
    genomegff = sys.argv[1]
    annotype = sys.argv[2]
    annokey = sys.argv[3]
    modified_gff = sys.argv[4]
    filtered_gff = sys.argv[5]
    overlapfilter(genomegff, annotype, annokey, modified_gff, filtered_gff)