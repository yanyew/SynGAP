#!/usr/bin/env python3
# -*- coding:utf-8 -*-
import os
import sys
import subprocess
import multiprocessing
from Bio import SeqIO


def miniprotcmd(miniprotdir, kmer1, kmer2, outs, intron, n):
    query = open('./query/' + str(n) + '.query.fa', 'r')
    subject = open('./subject/' + str(n) + '.subject.fa', 'r')
    blastlog = open('./miniprotoutput/' + str(n) + '.gapblast.log', 'a')
    queryseq = SeqIO.to_dict(SeqIO.parse(query, 'fasta'))
    subjectseq = SeqIO.to_dict(SeqIO.parse(subject, 'fasta'))
    len_list_query = []
    len_list_subject = []
    for key in queryseq:
        seq = queryseq[key]
        len_list_query.append(len(seq))
    for key in subjectseq:
        seq = subjectseq[key]
        len_list_subject.append(len(seq))
    if len_list_query and len_list_subject:
        if min(len_list_query) >= 7 and min(len_list_subject) >= 7:
            miniprotcommandline = str(miniprotdir) + '/miniprot -t1 --gff -k' + str(kmer1) + ' -l' + str(kmer2) + \
                                  ' --outs=' + str(outs) + ' -G' + str(intron) + \
                                  ' ./subject/' + str(n) + '.subject.fa' + \
                                  ' ./query/' + str(n) + '.query.fa > ./miniprotoutput/' + str(n) + '.gff'
            blastlog.write('command:' + str(miniprotcommandline) + '\n')
            p = subprocess.run(str(miniprotcommandline), shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                               encoding="utf-8")
            blastlog.write(p.stdout)
            blastlog.write('miniprot' + str(n) + 'finished' + '\n' + '\n')
    query.close()
    subject.close()
    blastlog.close()


def gapanno(miniprotdir, kmer1, kmer2, outs, intron, pnumber):
    inf_querynumber = open('./query/querynumber.txt', 'r')
    inf_subjectnumber = open('./subject/subjectnumber.txt', 'r')
    os.mkdir('./miniprotoutput')
    pool = multiprocessing.Pool(processes=int(pnumber))
    querynumber = inf_querynumber.readline()
    subjectnumber = inf_subjectnumber.readline()
    if querynumber != subjectnumber:
        print('Query number differences from subject number. Please check the previous step(gapseq).')
    else:
        for n in range(int(querynumber)):
            pool.apply_async(miniprotcmd, (miniprotdir, kmer1, kmer2, outs, intron, n))
    pool.close()
    pool.join()
    mergecommandline = 'find ./miniprotoutput -name "*.gff" | xargs -i cat {} | ' \
                       'sed "s/##gff-version 3//g" | sed "/^\s*$/d" > ' + '../modified.raw.gff'
    subprocess.run(str(mergecommandline), shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                   encoding="utf-8")
    inf_querynumber.close()
    inf_subjectnumber.close()


if __name__ == '__main__':
    miniprotdir = sys.argv[1]
    kmer1 = sys.argv[2]
    kmer2 = sys.argv[3]
    outs = sys.argv[4]
    intron = sys.argv[5]
    pnumber = sys.argv[6]
    gapanno(miniprotdir, kmer1, kmer2, outs, intron, pnumber)