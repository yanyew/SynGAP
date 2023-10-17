#!/usr/bin/env python3
# -*- coding:utf-8 -*-
import sys
import os

def reformatneedle(path, sp, n):
    seq_id_dict = {}
    Length_dict = {}
    Identity_dict = {}
    Similarity_dict = {}
    Gaps_dict = {}
    Score_dict = {}
    R_dict = {}
    for i in range(int(n)):
        nf = open(path + '/' + sp + '_global_alignment/' + str(i) + '.needle', 'r')
        seq_id1 = ''
        seq_id2 = ''
        Length = ''
        Identity_counts = ''
        Similarity_counts = ''
        Gaps_counts = ''
        Score = ''
        while True:
            line = nf.readline()
            if not line: break
            if line[0:4] == '# 1:':
                seq_id1 = line[5:].strip('\n')
            elif line[0:4] == '# 2:':
                seq_id2 = line[5:].strip('\n')
            elif line[0:9] == '# Length:':
                Length = line[10:].strip('\n')
            elif line[0:11] == '# Identity:':
                Identity_counts = ((line.split(':')[1]).split('/')[0]).replace(' ','')
            elif line[0:13] == '# Similarity:':
                Similarity_counts = ((line.split(':')[1]).split('/')[0]).replace(' ','')
            elif line[0:7] == '# Gaps:':
                Gaps_counts = ((line.split(':')[1]).split('/')[0]).replace(' ','')
            elif line[0:8] == '# Score:':
                Score = ((line.split(':')[1]).rstrip('\n')).replace(' ','')
        Identity = int(Identity_counts) / int(Length)
        Similarity = int(Similarity_counts) / int(Length)
        Gaps = int(Gaps_counts) / int(Length)
        R = Similarity * (1 - Gaps)
        seq_id_dict[seq_id1] = seq_id2
        Length_dict[seq_id1] = Length
        Identity_dict[seq_id1] = Identity
        Similarity_dict[seq_id1] = Similarity
        Gaps_dict[seq_id1] = Gaps
        Score_dict[seq_id1] = Score
        R_dict[seq_id1] = R
        nf.close()
    return(seq_id_dict, Length_dict, Identity_dict, Similarity_dict, Gaps_dict, Score_dict, R_dict)


def Rlabel(sp, gff, outR, outgff):
    path = os.getcwd()
    gf = open(gff, 'r')
    numf = open(path + '/' + sp + '_global_alignment/n.txt', 'r')
    oR = open(outR, 'w')
    og = open(outgff, 'w')
    oR.write("ID\torigin_ID\tLength\tIdentity\tSimilarity\tGaps\tScore\tR\tClassification\n")
    n = numf.readline()
    seq_id_dict, Length_dict, Identity_dict, Similarity_dict, Gaps_dict, Score_dict, R_dict = reformatneedle(path, sp, n)
    while True:
        line = gf.readline()
        if not line: break
        if line[0] != '#':
            type = line.split('\t')[2]
            if type == 'mRNA':
                attributes = line.split('\t')[8].rstrip()
                ID = (attributes.split(';')[0]).split('=')[1]
                if ID in seq_id_dict:
                    if R_dict[ID] >= 0.5:
                        oR.write(ID + '\t' + seq_id_dict[ID] + '\t' + Length_dict[ID] + '\t' +
                                 str(round(Identity_dict[ID], 4)) + '\t' + str(round(Similarity_dict[ID], 4)) + '\t' +
                                 str(round(Gaps_dict[ID], 4)) + '\t' + Score_dict[ID] + '\t' +
                                 str(round(R_dict[ID], 4)) + '\treal_gene\n')
                        og.write(line.rstrip('\n') + ';R=' + str(round(R_dict[ID], 4)) + ';classfication=real_gene\n')
                    elif R_dict[ID] < 0.5:
                        oR.write(ID + '\t' + seq_id_dict[ID] + '\t' + Length_dict[ID] + '\t' +
                                 str(round(Identity_dict[ID], 4)) + '\t' + str(round(Similarity_dict[ID], 4)) + '\t' +
                                 str(round(Gaps_dict[ID], 4)) + '\t' + Score_dict[ID] + '\t' +
                                 str(round(R_dict[ID], 4)) + '\tpotential_pseudogene\n')
                        og.write(line.rstrip('\n') + ';R=' + str(round(R_dict[ID], 4)) +
                                 ';classfication=potential_pseudogene\n')
                elif ID not in seq_id_dict:
                    og.write(line.rstrip('\n') + '\n')
            elif type == 'CDS':
                og.write(line.rstrip('\n') + '\n')
    numf.close()
    gf.close()
    oR.close()
    og.close()


if __name__ == '__main__':
    sp = sys.argv[1]
    gff = sys.argv[2]
    outR = sys.argv[3]
    outgff = sys.argv[4]
    Rlabel(sp, gff, outR, outgff)