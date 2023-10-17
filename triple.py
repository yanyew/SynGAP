#!/usr/bin/env python3
# -*- coding:utf-8 -*-
import os
import datetime
from scripts import SynGAP_genblastg, SynGAP_miniprot, redundantfilter


def triple(args):
    global workDir
    workDir = os.getcwd()
    mission_time = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    SynGAP_triple_Dir = workDir + '/SynGAP_triple_' + mission_time
    SynGAP_triple_results_Dir = SynGAP_triple_Dir + '/results'
    SynGAP_triple_round1Dir = SynGAP_triple_Dir + '/Round1'
    os.mkdir(SynGAP_triple_Dir)
    os.mkdir(SynGAP_triple_results_Dir)
    os.mkdir(SynGAP_triple_round1Dir)
    sp1fa_round1 = SynGAP_triple_round1Dir + '/' + str(args.sp1) + '.fa'
    sp1gff_round1 = SynGAP_triple_round1Dir + '/' + str(args.sp1) + '.gff3'
    sp2fa_round1 = SynGAP_triple_round1Dir + '/' + str(args.sp2) + '.fa'
    sp2gff_round1 = SynGAP_triple_round1Dir + '/' + str(args.sp2) + '.gff3'
    sp3fa_round1 = SynGAP_triple_round1Dir + '/' + str(args.sp3) + '.fa'
    sp3gff_round1 = SynGAP_triple_round1Dir + '/' + str(args.sp3) + '.gff3'

    # Start round 1
    print('[\033[0;36mINFO\033[0m] Renaming the input files, please wait ...')
    print('[\033[0;36mINFO\033[0m] '
          'Rename `\033[0;35m' + str(args.sp1fa) + '\033[0m` into `\033[0;35m' + sp1fa_round1 + '\033[0m`')
    os.system('cp ' + str(args.sp1fa) + ' ' + sp1fa_round1)
    print('[\033[0;36mINFO\033[0m] '
          'Rename `\033[0;35m' + str(args.sp1gff) + '\033[0m` into `\033[0;35m' + sp1gff_round1 + '\033[0m`')
    os.system('cp ' + str(args.sp1gff) + ' ' + sp1gff_round1)
    print('[\033[0;36mINFO\033[0m] '
          'Rename `\033[0;35m' + str(args.sp2fa) + '\033[0m` into `\033[0;35m' + sp2fa_round1 + '\033[0m`')
    os.system('cp ' + str(args.sp2fa) + ' ' + sp2fa_round1)
    print('[\033[0;36mINFO\033[0m] '
          'Rename `\033[0;35m' + str(args.sp2gff) + '\033[0m` into `\033[0;35m' + sp2gff_round1 + '\033[0m`')
    os.system('cp ' + str(args.sp2gff) + ' ' + sp2gff_round1)
    print('[\033[0;36mINFO\033[0m] '
          'Rename `\033[0;35m' + str(args.sp3fa) + '\033[0m` into `\033[0;35m' + sp3fa_round1 + '\033[0m`')
    os.system('cp ' + str(args.sp3fa) + ' ' + sp3fa_round1)
    print('[\033[0;36mINFO\033[0m] '
          'Rename `\033[0;35m' + str(args.sp3gff) + '\033[0m` into `\033[0;35m' + sp3gff_round1 + '\033[0m`')
    os.system('cp ' + str(args.sp3gff) + ' ' + sp3gff_round1)
    print('[\033[0;36mINFO\033[0m] Running Done!\n')
    os.chdir(SynGAP_triple_round1Dir)
    sp1_sp2_round1_resultDir = ''
    sp1_sp3_round1_resultDir = ''
    sp2_sp3_round1_resultDir = ''
    if args.process == 'genblastg':
        sp1_sp2_round1_resultDir = SynGAP_genblastg.SynGAP(
            args.sp1, args.sp2,
            args.annoType1, args.annoKey1, args.annoparentKey1,
            args.annoType2, args.annoKey2, args.annoparentKey2,
            args.datatype, args.cscore, args.threads, args.evalue, args.rank, args.coverage
        )
        sp1_sp3_round1_resultDir = SynGAP_genblastg.SynGAP(
            args.sp1, args.sp3,
            args.annoType1, args.annoKey1, args.annoparentKey1,
            args.annoType3, args.annoKey3, args.annoparentKey3,
            args.datatype, args.cscore, args.threads, args.evalue, args.rank, args.coverage
        )
        sp2_sp3_round1_resultDir = SynGAP_genblastg.SynGAP(
            args.sp2, args.sp3,
            args.annoType2, args.annoKey2, args.annoparentKey2,
            args.annoType3, args.annoKey3, args.annoparentKey3,
            args.datatype, args.cscore, args.threads, args.evalue, args.rank, args.coverage
        )
    elif args.process == 'miniprot':
        sp1_sp2_round1_resultDir = SynGAP_miniprot.SynGAP(
            args.sp1, args.sp2,
            args.annoType1, args.annoKey1, args.annoparentKey1,
            args.annoType2, args.annoKey2, args.annoparentKey2,
            args.datatype, args.cscore, args.threads, args.kmer1, args.kmer2, args.outs, args.intron
        )
        sp1_sp3_round1_resultDir = SynGAP_miniprot.SynGAP(
            args.sp1, args.sp3,
            args.annoType1, args.annoKey1, args.annoparentKey1,
            args.annoType3, args.annoKey3, args.annoparentKey3,
            args.datatype, args.cscore, args.threads, args.kmer1, args.kmer2, args.outs, args.intron
        )
        sp2_sp3_round1_resultDir = SynGAP_miniprot.SynGAP(
            args.sp2, args.sp3,
            args.annoType2, args.annoKey2, args.annoparentKey2,
            args.annoType3, args.annoKey3, args.annoparentKey3,
            args.datatype, args.cscore, args.threads, args.kmer1, args.kmer2, args.outs, args.intron
        )
    sp1modgff_c1_round1 = sp1_sp2_round1_resultDir + '/' + str(args.sp1) + '.modified.filtered.R.clusterd.gff'
    sp1modgff_c2_round1 = sp1_sp3_round1_resultDir + '/' + str(args.sp1) + '.modified.filtered.R.clusterd.gff'
    sp2modgff_c1_round1 = sp1_sp2_round1_resultDir + '/' + str(args.sp2) + '.modified.filtered.R.clusterd.gff'
    sp2modgff_c2_round1 = sp2_sp3_round1_resultDir + '/' + str(args.sp2) + '.modified.filtered.R.clusterd.gff'
    sp3modgff_c1_round1 = sp1_sp3_round1_resultDir + '/' + str(args.sp3) + '.modified.filtered.R.clusterd.gff'
    sp3modgff_c2_round1 = sp2_sp3_round1_resultDir + '/' + str(args.sp3) + '.modified.filtered.R.clusterd.gff'
    sp1modgff_c_round1 = SynGAP_triple_round1Dir + '/' + str(args.sp1) + '.SynGAPround1.c.gff'
    sp2modgff_c_round1 = SynGAP_triple_round1Dir + '/' + str(args.sp2) + '.SynGAPround1.c.gff'
    sp3modgff_c_round1 = SynGAP_triple_round1Dir + '/' + str(args.sp3) + '.SynGAPround1.c.gff'
    os.system('cat ' + sp1modgff_c1_round1 + ' ' + sp1modgff_c2_round1 + ' > ' + sp1modgff_c_round1)
    os.system('cat ' + sp2modgff_c1_round1 + ' ' + sp2modgff_c2_round1 + ' > ' + sp2modgff_c_round1)
    os.system('cat ' + sp3modgff_c1_round1 + ' ' + sp3modgff_c2_round1 + ' > ' + sp3modgff_c_round1)

    # Filter out the modified annotations that are redundant and overlapping with the annotations from genome .gff3
    sp1modbed_round1 = SynGAP_triple_round1Dir + '/filter/' + str(args.sp1) + '.SynGAPround1.c.bed'
    sp2modbed_round1 = SynGAP_triple_round1Dir + '/filter/' + str(args.sp2) + '.SynGAPround1.c.bed'
    sp3modbed_round1 = SynGAP_triple_round1Dir + '/filter/' + str(args.sp3) + '.SynGAPround1.c.bed'
    sp1modgff_disredundant_round1 = SynGAP_triple_round1Dir + '/filter/' + str(
        args.sp1) + '.SynGAPround1.I.disredundant.gff'
    sp2modgff_disredundant_round1 = SynGAP_triple_round1Dir + '/filter/' + str(
        args.sp2) + '.SynGAPround1.I.disredundant.gff'
    sp3modgff_disredundant_round1 = SynGAP_triple_round1Dir + '/filter/' + str(
        args.sp3) + '.SynGAPround1.I.disredundant.gff'
    sp1modgff_round1 = SynGAP_triple_round1Dir + '/' + str(args.sp1) + '.SynGAPround1.gff'
    sp2modgff_round1 = SynGAP_triple_round1Dir + '/' + str(args.sp2) + '.SynGAPround1.gff'
    sp3modgff_round1 = SynGAP_triple_round1Dir + '/' + str(args.sp3) + '.SynGAPround1.gff'
    print('[\033[0;36mINFO\033[0m] Filtering out the redundant modified annotations, please wait ...')
    os.mkdir(SynGAP_triple_round1Dir + '/filter')
    os.system(
        'python -m jcvi.formats.gff bed --type=mRNA --key=ID --parent_key=Parent ' + sp1modgff_c_round1 + ' -o ' + sp1modbed_round1)
    os.system(
        'python -m jcvi.formats.gff bed --type=mRNA --key=ID --parent_key=Parent ' + sp2modgff_c_round1 + ' -o ' + sp2modbed_round1)
    os.system(
        'python -m jcvi.formats.gff bed --type=mRNA --key=ID --parent_key=Parent ' + sp3modgff_c_round1 + ' -o ' + sp3modbed_round1)
    redundantfilter.redundantfilter(sp1modgff_c_round1, sp1modbed_round1, sp1modgff_disredundant_round1)
    print('[\033[0;36mINFO\033[0m] Running Done!')
    redundantfilter.redundantfilter(sp2modgff_c_round1, sp2modbed_round1, sp2modgff_disredundant_round1)
    print('[\033[0;36mINFO\033[0m] Running Done!')
    redundantfilter.redundantfilter(sp3modgff_c_round1, sp3modbed_round1, sp3modgff_disredundant_round1)
    print('[\033[0;36mINFO\033[0m] Running Done!')
    os.system('cat ' + sp1gff_round1 + ' ' + sp1modgff_disredundant_round1 + ' > ' + sp1modgff_round1)
    os.system('cat ' + sp2gff_round1 + ' ' + sp2modgff_disredundant_round1 + ' > ' + sp2modgff_round1)
    os.system('cat ' + sp3gff_round1 + ' ' + sp3modgff_disredundant_round1 + ' > ' + sp3modgff_round1)
    sp1modgff_SynGAPround1 = SynGAP_triple_results_Dir + '/' + str(args.sp1) + '.SynGAPround1.gff'
    sp2modgff_SynGAPround1 = SynGAP_triple_results_Dir + '/' + str(args.sp2) + '.SynGAPround1.gff'
    sp3modgff_SynGAPround1 = SynGAP_triple_results_Dir + '/' + str(args.sp3) + '.SynGAPround1.gff'
    os.system('ln -s ' + sp1modgff_round1 + ' ' + sp1modgff_SynGAPround1)
    os.system('ln -s ' + sp2modgff_round1 + ' ' + sp2modgff_SynGAPround1)
    os.system('ln -s ' + sp3modgff_round1 + ' ' + sp3modgff_SynGAPround1)

    # Start round 2
    os.chdir(workDir)
    SynGAP_triple_round2Dir = SynGAP_triple_Dir + '/Round2'
    os.mkdir(SynGAP_triple_round2Dir)
    sp1fa_round2 = SynGAP_triple_round2Dir + '/' + str(args.sp1) + '.fa'
    sp1gff_round2 = SynGAP_triple_round2Dir + '/' + str(args.sp1) + '.gff3'
    sp2fa_round2 = SynGAP_triple_round2Dir + '/' + str(args.sp2) + '.fa'
    sp2gff_round2 = SynGAP_triple_round2Dir + '/' + str(args.sp2) + '.gff3'
    sp3fa_round2 = SynGAP_triple_round2Dir + '/' + str(args.sp3) + '.fa'
    sp3gff_round2 = SynGAP_triple_round2Dir + '/' + str(args.sp3) + '.gff3'
    print('[\033[0;36mINFO\033[0m] Renaming the input files, please wait ...')
    print('[\033[0;36mINFO\033[0m] '
          'Rename `\033[0;35m' + str(args.sp1fa) + '\033[0m` into `\033[0;35m' + sp1fa_round2 + '\033[0m`')
    os.system('cp ' + str(args.sp1fa) + ' ' + sp1fa_round2)
    print('[\033[0;36mINFO\033[0m] '
          'Rename `\033[0;35m' + str(args.sp1gff) + '\033[0m` into `\033[0;35m' + sp1gff_round2 + '\033[0m`')
    os.system('cp ' + sp1modgff_SynGAPround1 + ' ' + sp1gff_round2)
    print('[\033[0;36mINFO\033[0m] '
          'Rename `\033[0;35m' + str(args.sp2fa) + '\033[0m` into `\033[0;35m' + sp2fa_round2 + '\033[0m`')
    os.system('cp ' + str(args.sp2fa) + ' ' + sp2fa_round2)
    print('[\033[0;36mINFO\033[0m] '
          'Rename `\033[0;35m' + str(args.sp2gff) + '\033[0m` into `\033[0;35m' + sp2gff_round2 + '\033[0m`')
    os.system('cp ' + sp2modgff_SynGAPround1 + ' ' + sp2gff_round2)
    print('[\033[0;36mINFO\033[0m] '
          'Rename `\033[0;35m' + str(args.sp3fa) + '\033[0m` into `\033[0;35m' + sp3fa_round2 + '\033[0m`')
    os.system('cp ' + str(args.sp3fa) + ' ' + sp3fa_round2)
    print('[\033[0;36mINFO\033[0m] '
          'Rename `\033[0;35m' + str(args.sp3gff) + '\033[0m` into `\033[0;35m' + sp3gff_round2 + '\033[0m`')
    os.system('cp ' + sp3modgff_SynGAPround1 + ' ' + sp3gff_round2)
    print('[\033[0;36mINFO\033[0m] Running Done!\n')
    os.chdir(SynGAP_triple_round2Dir)
    sp1_sp2_round2_resultDir = ''
    sp1_sp3_round2_resultDir = ''
    sp2_sp3_round2_resultDir = ''
    if args.process == 'genblastg':
        sp1_sp2_round2_resultDir = SynGAP_genblastg.SynGAP(
            args.sp1, args.sp2,
            args.annoType1, args.annoKey1, args.annoparentKey1,
            args.annoType2, args.annoKey2, args.annoparentKey2,
            args.datatype, args.cscore, args.threads, args.evalue, args.rank, args.coverage
        )
        sp1_sp3_round2_resultDir = SynGAP_genblastg.SynGAP(
            args.sp1, args.sp3,
            args.annoType1, args.annoKey1, args.annoparentKey1,
            args.annoType3, args.annoKey3, args.annoparentKey3,
            args.datatype, args.cscore, args.threads, args.evalue, args.rank, args.coverage
        )
        sp2_sp3_round2_resultDir = SynGAP_genblastg.SynGAP(
            args.sp2, args.sp3,
            args.annoType2, args.annoKey2, args.annoparentKey2,
            args.annoType3, args.annoKey3, args.annoparentKey3,
            args.datatype, args.cscore, args.threads, args.evalue, args.rank, args.coverage
        )
    elif args.process == 'miniprot':
        sp1_sp2_round2_resultDir = SynGAP_miniprot.SynGAP(
            args.sp1, args.sp2,
            args.annoType1, args.annoKey1, args.annoparentKey1,
            args.annoType2, args.annoKey2, args.annoparentKey2,
            args.datatype, args.cscore, args.threads, args.kmer1, args.kmer2, args.outs, args.intron
        )
        sp1_sp3_round2_resultDir = SynGAP_miniprot.SynGAP(
            args.sp1, args.sp3,
            args.annoType1, args.annoKey1, args.annoparentKey1,
            args.annoType3, args.annoKey3, args.annoparentKey3,
            args.datatype, args.cscore, args.threads, args.kmer1, args.kmer2, args.outs, args.intron
        )
        sp2_sp3_round2_resultDir = SynGAP_miniprot.SynGAP(
            args.sp2, args.sp3,
            args.annoType2, args.annoKey2, args.annoparentKey2,
            args.annoType3, args.annoKey3, args.annoparentKey3,
            args.datatype, args.cscore, args.threads, args.kmer1, args.kmer2, args.outs, args.intron
        )
    sp1modgff_I1_round2 = sp1_sp2_round2_resultDir + '/' + str(args.sp1) + '.modified.filtered.R.clusterd.gff'
    sp1modgff_I2_round2 = sp1_sp3_round2_resultDir + '/' + str(args.sp1) + '.modified.filtered.R.clusterd.gff'
    sp2modgff_c1_round2 = sp1_sp2_round2_resultDir + '/' + str(args.sp2) + '.modified.filtered.R.clusterd.gff'
    sp2modgff_c2_round2 = sp2_sp3_round2_resultDir + '/' + str(args.sp2) + '.modified.filtered.R.clusterd.gff'
    sp3modgff_c1_round2 = sp1_sp3_round2_resultDir + '/' + str(args.sp3) + '.modified.filtered.R.clusterd.gff'
    sp3modgff_c2_round2 = sp2_sp3_round2_resultDir + '/' + str(args.sp3) + '.modified.filtered.R.clusterd.gff'
    sp1modgff_c_round2 = SynGAP_triple_round2Dir + '/' + str(args.sp1) + '.SynGAPround2.c.gff'
    sp2modgff_c_round2 = SynGAP_triple_round2Dir + '/' + str(args.sp2) + '.SynGAPround2.c.gff'
    sp3modgff_c_round2 = SynGAP_triple_round2Dir + '/' + str(args.sp3) + '.SynGAPround2.c.gff'
    os.system('cat ' + sp1modgff_I1_round2 + ' ' + sp1modgff_I2_round2 + ' > ' + sp1modgff_c_round2)
    os.system('cat ' + sp2modgff_c1_round2 + ' ' + sp2modgff_c2_round2 + ' > ' + sp2modgff_c_round2)
    os.system('cat ' + sp3modgff_c1_round2 + ' ' + sp3modgff_c2_round2 + ' > ' + sp3modgff_c_round2)

    # Filter out the modified annotations that are redundant and overlapping with the annotations from genome .gff3
    sp1modbed_round2 = SynGAP_triple_round2Dir + '/filter/' + str(args.sp1) + '.SynGAPround2.c.bed'
    sp2modbed_round2 = SynGAP_triple_round2Dir + '/filter/' + str(args.sp2) + '.SynGAPround2.c.bed'
    sp3modbed_round2 = SynGAP_triple_round2Dir + '/filter/' + str(args.sp3) + '.SynGAPround2.c.bed'
    sp1modgff_disredundant_round2 = SynGAP_triple_round2Dir + '/filter/' + str(
        args.sp1) + '.SynGAPround2.I.disredundant.gff'
    sp2modgff_disredundant_round2 = SynGAP_triple_round2Dir + '/filter/' + str(
        args.sp2) + '.SynGAPround2.I.disredundant.gff'
    sp3modgff_disredundant_round2 = SynGAP_triple_round2Dir + '/filter/' + str(
        args.sp3) + '.SynGAPround2.I.disredundant.gff'
    sp1modgff_round2 = SynGAP_triple_round2Dir + '/' + str(args.sp1) + '.SynGAPround2.gff'
    sp2modgff_round2 = SynGAP_triple_round2Dir + '/' + str(args.sp2) + '.SynGAPround2.gff'
    sp3modgff_round2 = SynGAP_triple_round2Dir + '/' + str(args.sp3) + '.SynGAPround2.gff'
    print('[\033[0;36mINFO\033[0m] Filtering out the redundant modified annotations, please wait ...')
    os.mkdir(SynGAP_triple_round2Dir + '/filter')
    os.system(
        'python -m jcvi.formats.gff bed --type=mRNA --key=ID --parent_key=Parent ' + sp1modgff_c_round2 + ' -o ' + sp1modbed_round2)
    os.system(
        'python -m jcvi.formats.gff bed --type=mRNA --key=ID --parent_key=Parent ' + sp2modgff_c_round2 + ' -o ' + sp2modbed_round2)
    os.system(
        'python -m jcvi.formats.gff bed --type=mRNA --key=ID --parent_key=Parent ' + sp3modgff_c_round2 + ' -o ' + sp3modbed_round2)
    redundantfilter.redundantfilter(sp1modgff_c_round2, sp1modbed_round2, sp1modgff_disredundant_round2)
    print('[\033[0;36mINFO\033[0m] Running Done!')
    redundantfilter.redundantfilter(sp2modgff_c_round2, sp2modbed_round2, sp2modgff_disredundant_round2)
    print('[\033[0;36mINFO\033[0m] Running Done!')
    redundantfilter.redundantfilter(sp3modgff_c_round2, sp3modbed_round2, sp3modgff_disredundant_round2)
    print('[\033[0;36mINFO\033[0m] Running Done!')
    os.system('cat ' + sp1gff_round2 + ' ' + sp1modgff_disredundant_round2 + ' > ' + sp1modgff_round2)
    os.system('cat ' + sp2gff_round2 + ' ' + sp2modgff_disredundant_round2 + ' > ' + sp2modgff_round2)
    os.system('cat ' + sp3gff_round2 + ' ' + sp3modgff_disredundant_round2 + ' > ' + sp3modgff_round2)
    sp1modgff_SynGAPtriple = SynGAP_triple_results_Dir + '/' + str(args.sp1) + '.SynGAPtriple.gff'
    sp2modgff_SynGAPtriple = SynGAP_triple_results_Dir + '/' + str(args.sp2) + '.SynGAPtriple.gff'
    sp3modgff_SynGAPtriple = SynGAP_triple_results_Dir + '/' + str(args.sp3) + '.SynGAPtriple.gff'
    os.system('ln -s ' + sp1modgff_round2 + ' ' + sp1modgff_SynGAPtriple)
    os.system('ln -s ' + sp2modgff_round2 + ' ' + sp2modgff_SynGAPtriple)
    os.system('ln -s ' + sp3modgff_round2 + ' ' + sp3modgff_SynGAPtriple)

    print('[\033[0;36mINFO\033[0m] SynGAP analysis for \033[0;35m' + str(args.sp1) +
          '\033[0m, \033[0;35m' + str(args.sp2) + '\033[0m and \033[0;35m' + str(args.sp3) + '\033[0m Done!')
    print('[\033[0;36mINFO\033[0m] Please check the result files in `\033[0;35m' + SynGAP_triple_results_Dir + '\033[0m`\n')