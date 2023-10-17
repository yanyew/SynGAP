#!/usr/bin/env python3
# -*- coding:utf-8 -*-
import sys
import os
import shlex
import shutil
import subprocess
import multiprocessing

def modifiedalign2origin(sp, infasta1, infasta2):
    inf1 = open(infasta1, 'r')
    inf2 = open(infasta2, 'r')
    seqs1 = {}
    seqs2 = {}
    seq_id = ''
    seq_id1_dict = {}
    seq_id2_dict = {}
    path = os.getcwd()
    os.mkdir(sp + '_global_alignment')
    while True:
        line = inf1.readline()
        if not line: break
        if line[0] == '>':
            seq_id = line[1:].strip()
            seqs1[seq_id] = ''
        else:
            seqs1[seq_id] += line.strip()
    while True:
        line = inf2.readline()
        if not line: break
        if line[0] == '>':
            seq_id = line[1:].strip()
            seqs2[seq_id] = ''
        else:
            seqs2[seq_id] += line.strip()
    n = 0
    for seq_id2 in seqs2:
        for seq_id1 in seqs1:
            if seq_id2 in seq_id1:
                outf_query = open(path + '/' + sp + '_global_alignment/' + str(n) + '.query.fa', 'a')
                outf_subject = open(path + '/' + sp + '_global_alignment/' + str(n) + '.subject.fa', 'a')
                outf_query.write('>' + seq_id1 + '\n' + seqs1[seq_id1])
                outf_subject.write('>' + seq_id2 + '\n' + seqs2[seq_id2])
                seq_id1_dict[n] = seq_id1
                seq_id2_dict[n] = seq_id2
                outf_query.close()
                outf_subject.close()
                n += 1
    inf1.close()
    inf2.close()
    nfile = open(path + '/' + sp + '_global_alignment/n.txt', 'a')
    nfile.write(str(n))
    return n, seq_id1_dict, seq_id2_dict


def needlecmd(path, sp, n, seq_id1_dict, seq_id2_dict):
    alignlog = open(path + '/' + sp + '_global_alignment/' + str(n) + '.alignlog', 'a')
    needlecommandline = 'needle ' + path + '/' + sp + '_global_alignment/' + str(
        n) + '.query.fa ' + path + '/' + sp + '_global_alignment/' + str(
        n) + '.subject.fa -outfile ' + path + '/' + sp + '_global_alignment/' + str(
        n) + '.needle -gapopen 10.0 -gapextend 0.5 -sprotein1 Y -sprotein2 Y -sid1 ' + \
                        seq_id1_dict[n] + ' -sid2 ' + seq_id2_dict[n]
    alignlog.write('Command:' + needlecommandline + '\n')
    p = subprocess.run(shlex.split(needlecommandline), shell=False, stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT, encoding="utf-8")
    alignlog.write(p.stdout + '\n' + p.stderr)
    alignlog.close()


def needle(sp, infasta1, infasta2, pnumber):
    path = os.getcwd()
    n, seq_id1_dict, seq_id2_dict = modifiedalign2origin(sp, infasta1, infasta2)
    pool = multiprocessing.Pool(processes=int(pnumber))
    for i in range(int(n)):
        pool.apply_async(needlecmd, (path, sp, i, seq_id1_dict, seq_id2_dict))
    pool.close()
    pool.join()


if __name__ == '__main__':
    sp = sys.argv[1]
    infasta1 = sys.argv[2]
    infasta2 = sys.argv[3]
    pnumber = sys.argv[4]
    needle(sp, infasta1, infasta2, pnumber)