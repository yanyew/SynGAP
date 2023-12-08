#!/usr/bin/env python3
# -*- coding:utf-8 -*-
import sys


def map_genomegff(in_genomegff, sp, annotype, annoKey, annoparentKey):
    inf_genomegff = open(in_genomegff, 'r')
    n = 1
    gene_order = {}
    gene_repeat = {}
    gene_list = []
    gene_rename_map = {}
    mRNA_rename_map = {}
    while True:
        line = inf_genomegff.readline()
        if not line: break
        if line[0] != '#' and len(line.split('\t')) == 9:
            type = line.split('\t')[2]
            if type == annotype:
                attributes = (line.split('\t')[8]).rstrip('\n')
                attributes_dict = {}
                attributes_list = attributes.split(';')
                for a in attributes_list:
                    if a != '':
                        attributes_dict[a.split('=')[0]] = a.split('=')[1]
                if annoparentKey in attributes_dict:
                    gene = attributes_dict[annoparentKey]
                    if gene not in gene_list:
                        gene_order[gene] = n
                        gene_repeat[gene] = 1
                        n += 1
                        gene_list.append(gene)
                        gene_rename_map[gene] = sp + 'g' + str(gene_order[gene])
                        if annoKey in attributes_dict:
                            name = attributes_dict[annoKey]
                            mRNA_rename_map[name] = sp + 'g' + str(gene_order[gene]) + '.' + str(gene_repeat[gene])
                    elif gene in gene_list:
                        gene_repeat[gene] += 1
                        if annoKey in attributes_dict:
                            name = attributes_dict[annoKey]
                            mRNA_rename_map[name] = sp + 'g' + str(gene_order[gene]) + '.' + str(gene_repeat[gene])
                elif annoparentKey not in attributes_dict:
                    gene = attributes_dict[annoKey]
                    if gene not in gene_list:
                        gene_order[gene] = n
                        gene_repeat[gene] = 1
                        n += 1
                        gene_list.append(gene)
                        gene_rename_map[gene] = sp + 'g' + str(gene_order[gene])
                        if annoKey in attributes_dict:
                            name = attributes_dict[annoKey]
                            mRNA_rename_map[name] = sp + 'g' + str(gene_order[gene]) + '.' + str(gene_repeat[gene])
                    elif gene in gene_list:
                        gene_repeat[gene] += 1
                        if annoKey in attributes_dict:
                            name = attributes_dict[annoKey]
                            mRNA_rename_map[name] = sp + 'g' + str(gene_order[gene]) + '.' + str(gene_repeat[gene])
    inf_genomegff.close()
    return gene_rename_map, mRNA_rename_map


def rename_gff(in_genomegff, sp, annotype, annoKey, annoparentKey, out_genomegff, out_annotypemap):
    inf_gff = open(in_genomegff, 'r')
    outf_gff = open(out_genomegff, 'w')
    outf_annotypemap = open(out_annotypemap, 'w')
    print('[\033[0;36mINFO\033[0m] Load file `\033[0;35m' + in_genomegff + '\033[0m`')
    gene_rename_map, mRNA_rename_map = map_genomegff(in_genomegff, sp, annotype, annoKey, annoparentKey)
    print('[\033[0;36mINFO\033[0m] Renaming gene ids...')
    while True:
        line = inf_gff.readline()
        if not line: break
        if line[0] != '#' and len(line.split('\t')) == 9:
            type = line.split('\t')[2]
            if type == 'gene':
                attributes = (line.split('\t')[8]).rstrip('\n')
                attributes_dict = {}
                attributes_list = attributes.split(';')
                for a in attributes_list:
                    if a != '':
                        attributes_dict[a.split('=')[0]] = a.split('=')[1]
                if annoKey in attributes_dict:
                    ID = attributes_dict['ID']
                    if ID in gene_rename_map:
                        line_renamed = line.replace(ID, gene_rename_map[ID])
                        outf_gff.write(line_renamed)
            elif type == annotype:
                attributes = (line.split('\t')[8]).rstrip('\n')
                attributes_dict = {}
                attributes_list = attributes.split(';')
                for a in attributes_list:
                    if a != '':
                        attributes_dict[a.split('=')[0]] = a.split('=')[1]
                if annoparentKey in attributes_dict:
                    parent = attributes_dict[annoparentKey]
                    name = attributes_dict[annoKey]
                    line_renamed = line.replace(name, mRNA_rename_map[name])
                    line_renamed2 = line_renamed.replace(parent, gene_rename_map[parent])
                    outf_gff.write(line_renamed2)
                elif annoparentKey not in attributes_dict:
                    name = attributes_dict[annoKey]
                    line_renamed = line.replace(name, mRNA_rename_map[name])
                    outf_gff.write(line_renamed)
            else:
                attributes = (line.split('\t')[8]).rstrip('\n')
                attributes_dict = {}
                attributes_list = attributes.split(';')
                for a in attributes_list:
                    if a != '':
                        attributes_dict[a.split('=')[0]] = a.split('=')[1]
                if 'Parent' in attributes_dict:
                    parent = attributes_dict['Parent']
                    if parent in mRNA_rename_map:
                        line_renamed = line.replace(parent, mRNA_rename_map[parent])
                        outf_gff.write(line_renamed)
                    elif parent in gene_rename_map:
                        line_renamed = line.replace(parent, gene_rename_map[parent])
                        outf_gff.write(line_renamed)
    for mRNA in mRNA_rename_map:
        outf_annotypemap.write(mRNA + '\t' + mRNA_rename_map[mRNA] + '\n')
    inf_gff.close()
    outf_gff.close()
    outf_annotypemap.close()


if __name__ == '__main__':
    in_genomegff = sys.argv[1]
    sp = sys.argv[2]
    annotype = sys.argv[3]
    annoKey = sys.argv[4]
    annoparentKey = sys.argv[5]
    out_genomegff = sys.argv[6]
    out_annotypemap = sys.argv[7]
    rename_gff(in_genomegff, sp, annotype, annoKey, annoparentKey, out_genomegff, out_annotypemap)
