#!/usr/bin/env python3
# -*- coding:utf-8 -*-
import sys


def dict_bed(infile):
    inf_bed = open(infile, 'r')
    gene_chr = {}
    gene_start = {}
    gene_end = {}
    geneid_list = []
    while True:
        line = inf_bed.readline()
        if not line: break
        chr = line.split('\t')[0]
        start = line.split('\t')[1]
        end = line.split('\t')[2]
        gene_id = line.split('\t')[3]
        gene_chr[gene_id] = chr
        gene_start[gene_id] = start
        gene_end[gene_id] = end
        geneid_list.append(gene_id)
    inf_bed.close()
    return gene_chr, gene_start, gene_end, geneid_list


def dict_modgff(in_modgff):
    inf_modgff = open(in_modgff, 'r')
    gene_seqid = {}
    gene_start = {}
    gene_end = {}
    gene_cdslength = {}
    mRNAs = {}
    while True:
        line = inf_modgff.readline()
        if not line: break
        if line[0] != '#' and len(line.split('\t')) == 9:
            seqid = line.split('\t')[0]
            type = line.split('\t')[2]
            start = line.split('\t')[3]
            end = line.split('\t')[4]
            attributes = (line.split('\t')[8]).rstrip('\n')
            attributes_list = attributes.split(';')
            attributes_dict = {}
            for a in attributes_list:
                if a != '':
                    attributes_dict[a.split('=')[0]] = a.split('=')[1]
            if type == 'mRNA':
                name = attributes_dict['ID']
                if name not in mRNAs:
                    mRNAs[name] = ''
                    mRNAs[name] += line
                elif name in mRNAs:
                    mRNAs[name] += line
                gene_seqid[name] = seqid
                gene_start[name] = start
                gene_end[name] = end
            elif type == 'CDS':
                Parent = attributes_dict['Parent']
                if Parent in mRNAs:
                    mRNAs[Parent] += line
                elif Parent not in mRNAs:
                    mRNAs[Parent] = ''
                    mRNAs[Parent] += line
                if Parent not in gene_cdslength:
                    gene_cdslength[Parent] = 0
                    gene_cdslength[Parent] += int(end) - int(start)
                elif Parent in gene_cdslength:
                    gene_cdslength[Parent] += int(end) - int(start)
    inf_modgff.close()
    return gene_seqid, gene_start, gene_end, gene_cdslength, mRNAs


def redundantfilter(modified_gff, modified_bed, filtered_gff):
    print('[\033[0;36mINFO\033[0m] Load file `\033[0;35m' + modified_gff + '\033[0m`')
    mod_seqid, mod_start, mod_end, mod_cdslength, mRNAs = dict_modgff(modified_gff)
    print('[\033[0;36mINFO\033[0m] Load file `\033[0;35m' + modified_bed + '\033[0m`')
    gene_chr, gene_start, gene_end, geneid_list = dict_bed(modified_bed)
    outfgff = open(filtered_gff, 'w')
    n = 0
    block = []
    overlap_blocks = {}
    retain_genes = []
    print('[\033[0;36mINFO\033[0m] Filtering ...')
    for geneid in geneid_list:
        p = geneid_list.index(geneid)
        if p == 0:
            block.append(geneid)
            overlap_blocks[n] = block
        elif p > 0:
            geneid_p = geneid_list[p - 1]
            if gene_chr[geneid] == gene_chr[geneid_p]:
                if int(gene_start[geneid]) < int(gene_end[geneid_p]):
                    block.append(geneid)
                elif int(gene_start[geneid]) >= int(gene_end[geneid_p]):
                    block = []
                    block.append(geneid)
                    n += 1
                    overlap_blocks[n] = block
            elif gene_chr[geneid] != gene_chr[geneid_p]:
                block = []
                block.append(geneid)
                n += 1
                overlap_blocks[n] = block
    for block_number in overlap_blocks:
        b = overlap_blocks[block_number]
        if len(b) > 1:
            cdslen = []
            for geneid in b:
                cds_len = mod_cdslength[geneid]
                cdslen.append(cds_len)
            cdslen_max = max(cdslen)
            i = cdslen.index(cdslen_max)
            retain_genes.append(b[i])
        elif len(b) == 1:
            retain_genes.append(b[0])
    for geneid in retain_genes:
        outfgff.write(mRNAs[geneid])
    outfgff.close()


if __name__ == '__main__':
    modified_gff = sys.argv[1]
    modified_bed = sys.argv[2]
    filtered_gff = sys.argv[3]
    redundantfilter(modified_gff, modified_bed, filtered_gff)
