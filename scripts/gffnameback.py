#!/usr/bin/env python3
# -*- coding:utf-8 -*-
import sys


def dict_map(in_map):
    inf_map = open(in_map, 'r')
    map = {}
    map_re = {}
    while True:
        line = inf_map.readline()
        if not line: break
        if line != '\n':
            origin_id = line.split('\t')[0]
            renamed_id = (line.split('\t')[1]).rstrip('\n')
            map[renamed_id] = origin_id
            map_re[origin_id] = renamed_id
    inf_map.close()
    return map, map_re


def nameback_gff(in_modifiedgff, in_annotypemap, out_namebackgff):
    print('[\033[0;36mINFO\033[0m] Load file `\033[0;35m' + in_modifiedgff + '\033[0m`')
    inf_gff = open(in_modifiedgff, 'r')
    outf_gff = open(out_namebackgff, 'w')
    print('[\033[0;36mINFO\033[0m] Load file `' + in_annotypemap + '`')
    mRNAmap, mRNAmap_re = dict_map(in_annotypemap)
    print('[\033[0;36mINFO\033[0m] Naming back ...')
    while True:
        line = inf_gff.readline()
        if not line: break
        if line[0] != '#' and len(line.split('\t')) == 9:
            type = line.split('\t')[2]
            if type == 'mRNA':
                attributes = (line.split('\t')[8]).rstrip('\n')
                attributes_dict = {}
                attributes_list = attributes.split(';')
                for a in attributes_list:
                    if a != '':
                        attributes_dict[a.split('=')[0]] = a.split('=')[1]
                ID = attributes_dict['ID'].split('_')[1]
                line_nameback = line.replace(ID, mRNAmap[ID])
                outf_gff.write(line_nameback)
            elif type == 'CDS':
                attributes = (line.split('\t')[8]).rstrip('\n')
                attributes_dict = {}
                attributes_list = attributes.split(';')
                for a in attributes_list:
                    if a != '':
                        attributes_dict[a.split('=')[0]] = a.split('=')[1]
                parent = attributes_dict['Parent'].split('_')[1]
                line_nameback = line.replace(parent, mRNAmap[parent])
                outf_gff.write(line_nameback)
    inf_gff.close()
    outf_gff.close()


if __name__ == '__main__':
    in_modifiedgff = sys.argv[1]
    in_annotypemap = sys.argv[2]
    out_namebackgff = sys.argv[3]
    nameback_gff(in_modifiedgff, in_annotypemap, out_namebackgff)