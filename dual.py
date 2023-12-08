#!/usr/bin/env python
# -*- coding:utf-8 -*-
import datetime
import os

from scripts import SynGAP_genblastg, SynGAP_miniprot


def dual(args):
    global workDir
    workDir = os.getcwd()
    mission_time = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    SynGAP_dual_Dir = workDir + '/SynGAP_dual_' + mission_time
    os.mkdir(SynGAP_dual_Dir)
    sp1fa = SynGAP_dual_Dir + '/' + str(args.sp1) + '.fa'
    sp1gff = SynGAP_dual_Dir + '/' + str(args.sp1) + '.gff3'
    sp2fa = SynGAP_dual_Dir + '/' + str(args.sp2) + '.fa'
    sp2gff = SynGAP_dual_Dir + '/' + str(args.sp2) + '.gff3'
    print('[\033[0;36mINFO\033[0m] Renaming the input files, please wait ...')
    print('[\033[0;36mINFO\033[0m] '
          'Rename `\033[0;35m' + str(args.sp1fa) + '\033[0m` into `\033[0;35m' + sp1fa + '\033[0m`')
    os.system('cp ' + str(args.sp1fa) + ' ' + sp1fa)
    print('[\033[0;36mINFO\033[0m] '
          'Rename `\033[0;35m' + str(args.sp1gff) + '\033[0m` into `\033[0;35m' + sp1gff + '\033[0m`')
    os.system('cp ' + str(args.sp1gff) + ' ' + sp1gff)
    print('[\033[0;36mINFO\033[0m] '
          'Rename `\033[0;35m' + str(args.sp2fa) + '\033[0m` into `\033[0;35m' + sp2fa + '\033[0m`')
    os.system('cp ' + str(args.sp2fa) + ' ' + sp2fa)
    print('[\033[0;36mINFO\033[0m] '
          'Rename `\033[0;35m' + str(args.sp2gff) + '\033[0m` into `\033[0;35m' + sp2gff + '\033[0m`')
    os.system('cp ' + str(args.sp2gff) + ' ' + sp2gff)
    print('[\033[0;36mINFO\033[0m] Running Done!\n')
    os.chdir(SynGAP_dual_Dir)
    if args.process == 'genblastg':
        SynGAP_genblastg.SynGAP(
            args.sp1, args.sp2,
            args.annoType1, args.annoKey1, args.annoparentKey1,
            args.annoType2, args.annoKey2, args.annoparentKey2,
            args.datatype, args.cscore, args.threads, args.evalue, args.rank, args.coverage
        )
    elif args.process == 'miniprot':
        SynGAP_miniprot.SynGAP(
            args.sp1, args.sp2,
            args.annoType1, args.annoKey1, args.annoparentKey1,
            args.annoType2, args.annoKey2, args.annoparentKey2,
            args.datatype, args.cscore, args.threads, args.kmer1, args.kmer2, args.outs, args.intron
        )
