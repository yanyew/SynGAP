#!/usr/bin/env python3
# -*- coding:utf-8 -*-
import sys


def dict_map(in_map):
    inf_map = open(in_map, 'r')
    map = {}
    while True:
        line = inf_map.readline()
        if not line: break
        if line.rstrip('\n') != '':
            map[line.rstrip('\n').split('\t')[0]] = line.rstrip('\n').split('\t')[1]
    inf_map.close()
    return map


def rename_anchors(in_anchors, sp1mRNAmap, sp2mRNAmap, out_anchors):
    inf_anchors = open(in_anchors, 'r')
    outf_anchors = open(out_anchors, 'w')
    print('[\033[0;36mINFO\033[0m] Renaming gene ids...')
    sp1map = dict_map(sp1mRNAmap)
    sp2map = dict_map(sp2mRNAmap)
    while True:
        line = inf_anchors.readline()
        if not line: break
        if line[0] == '#':
            outf_anchors.write(line)
        if line[0] != '#' and line.rstrip('\n') != '':
            sp1id = line.split('\t')[0]
            sp2id = line.split('\t')[1]
            outf_anchors.write(sp1map[sp1id] + '\t' + sp2map[sp2id] + '\n')
    inf_anchors.close()
    outf_anchors.close()


if __name__ == '__main__':
    in_anchors = sys.argv[1]
    sp1mRNAmap = sys.argv[2]
    sp2mRNAmap = sys.argv[3]
    out_anchors = sys.argv[4]
    rename_anchors(in_anchors, sp1mRNAmap, sp2mRNAmap, out_anchors)