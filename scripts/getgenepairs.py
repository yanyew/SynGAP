#!/usr/bin/env python3
# -*- coding:utf-8 -*-
import sys


def get_genepair(infile, outfile):
    inf = open(infile, 'r')
    outf = open(outfile, 'w')
    lines = inf.readlines()
    for line in lines:
        if line[0] != '#':
            outf.write(line.split('\t')[0] + '\t' + line.split('\t')[1] + '\n')
    inf.close()
    outf.close()


if __name__ == '__main__':
    infile = sys.argv[1]
    outfile = sys.argv[2]
    get_genepair(infile, outfile)