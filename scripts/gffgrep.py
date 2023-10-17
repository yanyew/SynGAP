#!/usr/bin/env python3
# -*- coding:utf-8 -*-
import sys


def gffgrep(gff, idI, idIIa, idIIb, idIII, outgffI, outgffIIa, outgffIIb, outgffIII):
    inf_gff = open(gff, 'r')
    inf_idI = open(idI, 'r')
    inf_idIIa = open(idIIa, 'r')
    inf_idIIb = open(idIIb, 'r')
    inf_idIII = open(idIII, 'r')
    outf_gffI = open(outgffI, 'w')
    outf_gffIIa = open(outgffIIa, 'w')
    outf_gffIIb = open(outgffIIb, 'w')
    outf_gffIII = open(outgffIII, 'w')
    idI_list = []
    idIIa_list = []
    idIIb_list = []
    idIII_list = []
    idI_lines = inf_idI.readlines()
    for line in idI_lines:
        idI_list.append(line.rstrip('\n'))
    idIIa_lines = inf_idIIa.readlines()
    for line in idIIa_lines:
        idIIa_list.append(line.rstrip('\n'))
    idIIb_lines = inf_idIIb.readlines()
    for line in idIIb_lines:
        idIIb_list.append(line.rstrip('\n'))
    idIII_lines = inf_idIII.readlines()
    for line in idIII_lines:
        idIII_list.append(line.rstrip('\n'))
    while True:
        line = inf_gff.readline()
        if not line: break
        if line[0] != '#' and line != '':
            type = line.split('\t')[2]
            attributes = (line.split('\t')[8]).rstrip('\n')
            attributes_list = attributes.split(';')
            attributes_dict = {}
            for a in attributes_list:
                if a != '':
                    attributes_dict[a.split('=')[0]] = a.split('=')[1]
            if type == 'mRNA':
                Target = attributes_dict['Target'].split(' ')[0]
                if Target in idI_list:
                    outf_gffI.write(line)
                elif Target in idIIa_list:
                    outf_gffIIa.write(line)
                elif Target in idIIb_list:
                    outf_gffIIb.write(line)
                if Target in idIII_list:
                    outf_gffIII.write(line)
            elif type == 'CDS':
                Target = attributes_dict['Target'].split(' ')[0]
                if Target in idI_list:
                    outf_gffI.write(line)
                elif Target in idIIa_list:
                    outf_gffIIa.write(line)
                elif Target in idIIb_list:
                    outf_gffIIb.write(line)
                elif Target in idIII_list:
                    outf_gffIII.write(line)
    inf_gff.close()
    inf_idI.close()
    inf_idIIa.close()
    inf_idIIb.close()
    inf_idIII.close()
    outf_gffI.close()
    outf_gffIIa.close()
    outf_gffIIb.close()
    outf_gffIII.close()


if __name__ == '__main__':
    gff = sys.argv[1]
    idI = sys.argv[2]
    idIIa = sys.argv[3]
    idIIb = sys.argv[4]
    idIII = sys.argv[5]
    outgffI = sys.argv[6]
    outgffIIa = sys.argv[7]
    outgffIIb = sys.argv[8]
    outgffIII = sys.argv[9]
    gffgrep(gff, idI, idIIa, idIIb, idIII, outgffI, outgffIIa, outgffIIb, outgffIII)