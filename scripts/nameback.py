#!/usr/bin/env python3
# -*- coding:utf-8 -*-
import sys


def dict_map(in_map):
    inf_map = open(in_map, 'r')
    map = {}
    while True:
        line = inf_map.readline()
        if not line: break
        if line != '\n':
            origin_id = line.split('\t')[0]
            renamed_id = (line.split('\t')[1]).rstrip('\n')
            map[renamed_id] = origin_id
    inf_map.close()
    return map


def nameback(infile, sp1mRNAmap, sp2mRNAmap, outfile):
    inf = open(infile, 'r')
    outf = open(outfile, 'w')
    print('[\033[0;36mINFO\033[0m] Load file `\033[0;35m' + infile + '\033[0m`')
    print('[\033[0;36mINFO\033[0m] Load file `' + sp1mRNAmap + '`')
    sp1map = dict_map(sp1mRNAmap)
    print('[\033[0;36mINFO\033[0m] Load file `' + sp2mRNAmap + '`')
    sp2map = dict_map(sp2mRNAmap)
    print('[\033[0;36mINFO\033[0m] Naming back ...')
    while True:
        line = inf.readline()
        if not line: break
        if line == '###\n' or line == '#\n':
            outf.write(line)
        elif line != '###\n' and line != '#\n':
            line_list = (line.rstrip('\n')).split('\t')
            if len(line_list) == 4:
                line_nameback = line.replace(line_list[0], sp1map[line_list[0]])
                line_nameback = line_nameback.replace(line_list[1], sp1map[line_list[1]])
                line_nameback = line_nameback.replace(line_list[2], sp2map[line_list[2]])
                line_nameback = line_nameback.replace(line_list[3], sp2map[line_list[3]])
                outf.write(line_nameback)
            if len(line_list) == 3:
                line_nameback = line.replace(line_list[0], sp1map[line_list[0]])
                line_nameback = line_nameback.replace(line_list[1], sp2map[line_list[1]])
                outf.write(line_nameback)
    inf.close()
    outf.close()


if __name__ == '__main__':
    infile = sys.argv[1]
    sp1mRNAmap = sys.argv[2]
    sp2mRNAmap = sys.argv[3]
    outfile = sys.argv[4]
    nameback(infile, sp1mRNAmap, sp2mRNAmap, outfile)
