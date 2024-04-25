#!/usr/bin/env python3
# -*- coding:utf-8 -*-
import os
import sys


def dict_bed(infile):
    inf_bed = open(infile, 'r')
    bed = {}
    while True:
        line = inf_bed.readline()
        if not line: break
        chr = line.split('\t')[0]
        start = line.split('\t')[1]
        end = line.split('\t')[2]
        id = line.split('\t')[3]
        if id not in bed:
            bed[id] = (chr, int(start), int(end))
    inf_bed.close()
    return bed


def dict_blockRbed(infile):
    inf_blockRbed = open(infile, 'r')
    blockR_bed = {}
    while True:
        line = inf_blockRbed.readline()
        if not line: break
        chr = line.split('\t')[0]
        start = line.split('\t')[1]
        end = line.split('\t')[2]
        blockR_id = line.split('\t')[3]
        R = line.split('\t')[4].rstrip('\n')
        if blockR_id not in blockR_bed:
            blockR_bed[blockR_id] = (chr, int(start), int(end), float(R))
    inf_blockRbed.close()
    return blockR_bed


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
                Identity_counts = ((line.split(':')[1]).split('/')[0]).replace(' ', '')
            elif line[0:13] == '# Similarity:':
                Similarity_counts = ((line.split(':')[1]).split('/')[0]).replace(' ', '')
            elif line[0:7] == '# Gaps:':
                Gaps_counts = ((line.split(':')[1]).split('/')[0]).replace(' ', '')
            elif line[0:8] == '# Score:':
                Score = ((line.split(':')[1]).rstrip('\n')).replace(' ', '')
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
    return (seq_id_dict, Length_dict, Identity_dict, Similarity_dict, Gaps_dict, Score_dict, R_dict)


def Rlabel(sp, spblockRbed, spnbed, spnblockRbed, gff, outR, outgff):
    path = os.getcwd()
    blockR_bed = dict_blockRbed(spblockRbed)
    nblockR_bed = dict_blockRbed(spnblockRbed)
    spn_bed = dict_bed(spnbed)
    gf = open(gff, 'r')
    numf = open(path + '/' + sp + '_global_alignment/n.txt', 'r')
    oR = open(outR, 'w')
    og = open(outgff, 'w')
    oglowR = open(gff.replace('.gff', '.lowR.gff'), 'w')
    oR.write("ID\torigin_ID\tLength\tIdentity\tSimilarity\tGaps\tScore\tR\tClassification\n")
    n = numf.readline()
    seq_id_dict, Length_dict, Identity_dict, Similarity_dict, Gaps_dict, Score_dict, R_dict = \
        reformatneedle(path, sp, n)
    lowR_list = []
    ID_block = {}
    Target_block = {}

    while True:
        line = gf.readline()
        if not line: break
        if line[0] != '#':
            type = line.split('\t')[2]
            if type == 'mRNA':
                chr = line.split('\t')[0].rstrip()
                start = int(line.split('\t')[3].rstrip())
                end = int(line.split('\t')[4].rstrip())
                attributes = (line.split('\t')[8]).rstrip('\n')
                attributes_list = attributes.split(';')
                attributes_dict = {}
                for a in attributes_list:
                    if a != '':
                        attributes_dict[a.split('=')[0]] = a.split('=')[1]
                ID = attributes_dict['ID']
                Target = attributes_dict['Target'].split(' ')[0]
                chr_Target = spn_bed[Target][0]
                start_Target = spn_bed[Target][1] + 1
                end_Target = spn_bed[Target][2]
                ID_block[ID] = {}
                Target_block[Target] = {}
                if ID in seq_id_dict:
                    for block in blockR_bed:
                        if chr == blockR_bed[block][0]:
                            if start <= blockR_bed[block][2] and end >= blockR_bed[block][1] + 1:
                                blockR = blockR_bed[block][3]
                                ID_block[ID][block] = blockR
                    for block in nblockR_bed:
                        if chr_Target == nblockR_bed[block][0]:
                            if start_Target <= nblockR_bed[block][2] and end_Target >= nblockR_bed[block][1] + 1:
                                blockR = nblockR_bed[block][3]
                                Target_block[Target][block] = blockR
                    block = ''
                    for i in ID_block[ID].keys():
                        if i in Target_block[Target].keys():
                            block = i
                    if block != '':
                        blockR = blockR_bed[block][3]
                        if R_dict[ID] >= blockR:
                            oR.write(ID + '\t' + seq_id_dict[ID] + '\t' + Length_dict[ID] + '\t' +
                                     str(round(Identity_dict[ID], 4)) + '\t' + str(round(Similarity_dict[ID], 4)) + '\t' +
                                     str(round(Gaps_dict[ID], 4)) + '\t' + Score_dict[ID] + '\t' +
                                     str(round(R_dict[ID], 4)) + '\treal_gene\n')
                            og.write(line.rstrip('\n') + ';R=' + str(round(R_dict[ID], 4)) + ';classfication=real_gene\n')
                        elif R_dict[ID] < blockR:
                            oR.write(ID + '\t' + seq_id_dict[ID] + '\t' + Length_dict[ID] + '\t' +
                                     str(round(Identity_dict[ID], 4)) + '\t' + str(round(Similarity_dict[ID], 4)) + '\t' +
                                     str(round(Gaps_dict[ID], 4)) + '\t' + Score_dict[ID] + '\t' +
                                     str(round(R_dict[ID], 4)) + '\tpotential_pseudogene\n')
                            oglowR.write(line.rstrip('\n') + ';R=' + str(round(R_dict[ID], 4)) + ';classfication=potential_pseudogene\n')
                            lowR_list.append(ID)
                    elif block == '':
                        if len(list(ID_block[ID].keys())) != 0:
                            block = list(ID_block[ID].keys())[0]
                            blockR = blockR_bed[block][3]
                            if R_dict[ID] >= blockR:
                                oR.write(ID + '\t' + seq_id_dict[ID] + '\t' + Length_dict[ID] + '\t' +
                                         str(round(Identity_dict[ID], 4)) + '\t' + str(
                                    round(Similarity_dict[ID], 4)) + '\t' +
                                         str(round(Gaps_dict[ID], 4)) + '\t' + Score_dict[ID] + '\t' +
                                         str(round(R_dict[ID], 4)) + '\treal_gene\n')
                                og.write(
                                    line.rstrip('\n') + ';R=' + str(round(R_dict[ID], 4)) + ';classfication=real_gene\n')
                            elif R_dict[ID] < blockR:
                                oR.write(ID + '\t' + seq_id_dict[ID] + '\t' + Length_dict[ID] + '\t' +
                                         str(round(Identity_dict[ID], 4)) + '\t' + str(
                                    round(Similarity_dict[ID], 4)) + '\t' +
                                         str(round(Gaps_dict[ID], 4)) + '\t' + Score_dict[ID] + '\t' +
                                         str(round(R_dict[ID], 4)) + '\tpotential_pseudogene\n')
                                oglowR.write(line.rstrip('\n') + ';R=' + str(
                                    round(R_dict[ID], 4)) + ';classfication=potential_pseudogene\n')
                                lowR_list.append(ID)
                        elif len(list(ID_block[ID].keys())) == 0 and len(list(Target_block[Target].keys())) != 0:
                            block = list(Target_block[Target].keys())[0]
                            blockR = blockR_bed[block][3]
                            if R_dict[ID] >= blockR:
                                oR.write(ID + '\t' + seq_id_dict[ID] + '\t' + Length_dict[ID] + '\t' +
                                         str(round(Identity_dict[ID], 4)) + '\t' + str(
                                    round(Similarity_dict[ID], 4)) + '\t' +
                                         str(round(Gaps_dict[ID], 4)) + '\t' + Score_dict[ID] + '\t' +
                                         str(round(R_dict[ID], 4)) + '\treal_gene\n')
                                og.write(
                                    line.rstrip('\n') + ';R=' + str(round(R_dict[ID], 4)) + ';classfication=real_gene\n')
                            elif R_dict[ID] < blockR:
                                oR.write(ID + '\t' + seq_id_dict[ID] + '\t' + Length_dict[ID] + '\t' +
                                         str(round(Identity_dict[ID], 4)) + '\t' + str(
                                    round(Similarity_dict[ID], 4)) + '\t' +
                                         str(round(Gaps_dict[ID], 4)) + '\t' + Score_dict[ID] + '\t' +
                                         str(round(R_dict[ID], 4)) + '\tpotential_pseudogene\n')
                                oglowR.write(line.rstrip('\n') + ';R=' + str(
                                    round(R_dict[ID], 4)) + ';classfication=potential_pseudogene\n')
                                lowR_list.append(ID)
                        elif len(list(ID_block[ID].keys())) == 0 and len(list(Target_block[Target].keys())) == 0:
                            blockR = 0.5
                            if R_dict[ID] >= blockR:
                                oR.write(ID + '\t' + seq_id_dict[ID] + '\t' + Length_dict[ID] + '\t' +
                                         str(round(Identity_dict[ID], 4)) + '\t' + str(
                                    round(Similarity_dict[ID], 4)) + '\t' +
                                         str(round(Gaps_dict[ID], 4)) + '\t' + Score_dict[ID] + '\t' +
                                         str(round(R_dict[ID], 4)) + '\treal_gene\n')
                                og.write(
                                    line.rstrip('\n') + ';R=' + str(round(R_dict[ID], 4)) + ';classfication=real_gene\n')
                            elif R_dict[ID] < blockR:
                                oR.write(ID + '\t' + seq_id_dict[ID] + '\t' + Length_dict[ID] + '\t' +
                                         str(round(Identity_dict[ID], 4)) + '\t' + str(
                                    round(Similarity_dict[ID], 4)) + '\t' +
                                         str(round(Gaps_dict[ID], 4)) + '\t' + Score_dict[ID] + '\t' +
                                         str(round(R_dict[ID], 4)) + '\tpotential_pseudogene\n')
                                oglowR.write(line.rstrip('\n') + ';R=' + str(
                                    round(R_dict[ID], 4)) + ';classfication=potential_pseudogene\n')
                                lowR_list.append(ID)
                elif ID not in seq_id_dict:
                    og.write(line.rstrip('\n') + '\n')
            elif type == 'CDS':
                attributes = (line.split('\t')[8]).rstrip('\n')
                attributes_list = attributes.split(';')
                attributes_dict = {}
                for a in attributes_list:
                    if a != '':
                        attributes_dict[a.split('=')[0]] = a.split('=')[1]
                if attributes_dict['Parent'] in lowR_list:
                    oglowR.write(line.rstrip('\n') + '\n')
                else:
                    og.write(line.rstrip('\n') + '\n')
    numf.close()
    gf.close()
    oR.close()
    og.close()
    oglowR.close()


if __name__ == '__main__':
    sp = sys.argv[1]
    spblockRbed = sys.argv[2]
    spnbed = sys.argv[3]
    spnblockRbed = sys.argv[4]
    gff = sys.argv[5]
    outR = sys.argv[6]
    outgff = sys.argv[7]
    Rlabel(sp, spblockRbed, spnbed, spnblockRbed, gff, outR, outgff)
