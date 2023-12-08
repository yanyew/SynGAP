#!/usr/bin/env python
# -*- coding:utf-8 -*-
import datetime
import os
import sys
from scripts import SynGAP_genblastg, SynGAP_miniprot, synteny


def get_synteny_number(file_path):
    column_1_list = []
    with open(file_path, 'r') as file:
        for line in file:
            if line[0] != '#':
                columns = line.strip().split('\t')
                column_1_list.append(columns[0])
    unique_elements = list(set(column_1_list))
    element_count = len(unique_elements)
    return unique_elements, element_count


def find_max_value_key(dictionary):
    max_value = max(dictionary.values())
    for key, value in dictionary.items():
        if value == max_value:
            return key


def master(args):
    script_dir, filename = os.path.split(os.path.abspath(sys.argv[0]))
    masterdb_dir = str(script_dir) + '/bin/masterdb/' + str(args.sp)
    splist = masterdb_dir + '/species.list'
    global sp_list
    if os.path.exists(splist):
        with open(splist, 'r') as file:
            lines = file.readlines()
        sp_list = [line.split('\t')[0] for line in lines]
    else:
        # print("The masterdb is incomplete. Please run `\033[0;35m'syngap init'\033[0m` first.")
        print("The masterdb is incomplete.\n"
              "Please check the manual in https://github.com/yanyew/SynGAP to get the download link.\n"
              "And then import the downloaded .tar.gz by running `\033[0;35m'syngap initdb'\033[0m`.")
        sys.exit()
    for sp in sp_list:
        if os.path.exists(masterdb_dir + '/' + sp + '.fa') and os.path.exists(masterdb_dir + '/' + sp + '.gff3'):
            continue
        else:
            # print("The masterdb is incomplete. Please run `\033[0;35m'syngap init'\033[0m` first.")
            print("The masterdb is incomplete.\n"
                  "Please check the manual in https://github.com/yanyew/SynGAP to get the download link.\n"
                  "And then import the downloaded .tar.gz by running `\033[0;35m'syngap initdb'\033[0m`.")
            sys.exit()
    print('[\033[0;36mINFO\033[0m] Check masterbd done! '
          'You can check the species list in \033[0;35m' + splist + '\033[0m.')

    global workDir
    workDir = os.getcwd()
    mission_time = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    SynGAP_master_Dir = workDir + '/SynGAP_master_' + mission_time
    SynGAP_synteny_Dir = SynGAP_master_Dir + '/optimal_synteny_search'
    os.mkdir(SynGAP_master_Dir)
    os.mkdir(SynGAP_synteny_Dir)
    os.chdir(SynGAP_synteny_Dir)
    print('[\033[0;36mINFO\033[0m] Searching the species keep optimal synteny with \033[0;35m' + args.sp1 +
          '\033[0m in masterdb, please wait ...')
    synteny_count = {}
    for sp in sp_list:
        sp_synteny_Dir = SynGAP_synteny_Dir + '/' + str(args.sp1) + '_' + sp
        os.mkdir(sp_synteny_Dir)
        os.chdir(sp_synteny_Dir)
        spfa = masterdb_dir + '/' + sp + '.fa'
        spgff = masterdb_dir + '/' + sp + '.gff3'
        sp1fa = sp_synteny_Dir + '/' + str(args.sp1) + '.fa'
        sp1gff = sp_synteny_Dir + '/' + str(args.sp1) + '.gff3'
        sp2fa = sp_synteny_Dir + '/' + sp + '.fa'
        sp2gff = sp_synteny_Dir + '/' + sp + '.gff3'
        os.symlink(workDir + '/' + str(args.sp1fa), sp1fa)
        os.symlink(workDir + '/' + str(args.sp1gff), sp1gff)
        os.symlink(spfa, sp2fa)
        os.symlink(spgff, sp2gff)
        anchors = synteny.synteny(
            args.sp1, sp,
            args.annoType1, args.annoKey1, args.annoparentKey1,
            'mRNA', 'ID', 'Parent',
            args.datatype, args.cscore, args.threads)
        unique_elements, element_count = get_synteny_number(anchors)
        synteny_count[sp] = element_count
    synteny_count_file = SynGAP_synteny_Dir + '/synteny_count.xls'
    with open(synteny_count_file, 'w') as file:
        for key, value in synteny_count.items():
            line = f"{key}\t{value}\n"
            file.write(line)
    sp2 = find_max_value_key(synteny_count)
    print('[\033[0;36mINFO\033[0m] Running Done!\n')

    os.chdir(SynGAP_master_Dir)
    spfa = masterdb_dir + '/' + sp2 + '.fa'
    spgff = masterdb_dir + '/' + sp2 + '.gff3'
    sp1fa = SynGAP_master_Dir + '/' + str(args.sp1) + '.fa'
    sp1gff = SynGAP_master_Dir + '/' + str(args.sp1) + '.gff3'
    sp2fa = SynGAP_master_Dir + '/' + sp2 + '.fa'
    sp2gff = SynGAP_master_Dir + '/' + sp2 + '.gff3'
    os.symlink(workDir + '/' + str(args.sp1fa), sp1fa)
    os.symlink(workDir + '/' + str(args.sp1gff), sp1gff)
    os.symlink(spfa, sp2fa)
    os.symlink(spgff, sp2gff)
    global SynGAP_results_Dir
    if args.process == 'genblastg':
        SynGAP_results_Dir = SynGAP_genblastg.SynGAP(
            args.sp1, sp2,
            args.annoType1, args.annoKey1, args.annoparentKey1,
            'mRNA', 'ID', 'Parent',
            args.datatype, args.cscore, args.threads, args.evalue, args.rank, args.coverage
        )
    elif args.process == 'miniprot':
        SynGAP_results_Dir = SynGAP_miniprot.SynGAP(
            args.sp1, sp2,
            args.annoType1, args.annoKey1, args.annoparentKey1,
            'mRNA', 'ID', 'Parent',
            args.datatype, args.cscore, args.threads, args.kmer1, args.kmer2, args.outs, args.intron
        )