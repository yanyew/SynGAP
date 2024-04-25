#!/usr/bin/env python
# -*- coding:utf-8 -*-
import multiprocessing
import os
import subprocess
import sys
import glob
from Bio import SeqIO


def genblastgcmd(genblastdir, evalue, rank, coverage, n):
    query = open('./query/' + str(n) + '.query.fa', 'r')
    subject = open('./subject/' + str(n) + '.subject.fa', 'r')
    blastlog = open('./genblastgoutput/' + str(n) + '.gapblast.log', 'a')
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
            genblastgcommandline = str(genblastdir) + '/genblast_v138_linux_x86_64 -p genblastg -P blast -q ./query/' \
                                   + str(n) + '.query.fa' + \
                                   ' -t ./subject/' + str(n) + '.subject.fa' + ' -e ' + str(evalue) + \
                                   ' -r ' + str(rank) + ' -g T -f F -c ' + str(coverage) + \
                                   ' -norepair -gff -o ./genblastgoutput/' + str(n) + '.gblast'
            blastlog.write('command:' + str(genblastgcommandline) + '\n')
            p = subprocess.run(str(genblastgcommandline), shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                               encoding="utf-8")
            blastlog.write(p.stdout)
            blastlog.write('genblastg' + str(n) + 'finished' + '\n' + '\n')
    query.close()
    subject.close()
    blastlog.close()


def gapanno(genblastdir, evalue, rank, coverage, pnumber):
    inf_querynumber = open('./query/querynumber.txt', 'r')
    inf_subjectnumber = open('./subject/subjectnumber.txt', 'r')
    os.mkdir('./genblastgoutput')
    pool = multiprocessing.Pool(processes=int(pnumber))
    querynumber = inf_querynumber.readline()
    subjectnumber = inf_subjectnumber.readline()
    if querynumber != subjectnumber:
        print('Query number is different from subject number. Please check the previous step (gapseq).')
    else:
        for n in range(int(querynumber)):
            pool.apply_async(genblastgcmd, (genblastdir, evalue, rank, coverage, n))
    pool.close()
    pool.join()
    directory = "./genblastgoutput"
    file_extension = "*.gff"
    file_list = glob.glob(os.path.join(directory, file_extension))
    merged_content = ""
    for file_path in file_list:
        with open(file_path, 'r') as file:
            file_content = file.read()
            merged_content += file_content
    modified_content = merged_content.replace("##gff-version 3\n", "").strip('')
    output_path = "../modified.raw.gff"
    with open(output_path, 'w') as output_file:
        output_file.write(modified_content)
    inf_querynumber.close()
    inf_subjectnumber.close()


if __name__ == '__main__':
    genblastdir = sys.argv[1]
    evalue = sys.argv[2]
    rank = sys.argv[3]
    coverage = sys.argv[4]
    pnumber = sys.argv[5]
    gapanno(genblastdir, evalue, rank, coverage, pnumber)
