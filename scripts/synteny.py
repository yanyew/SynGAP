#!/usr/bin/env python3
# -*- coding:utf-8 -*-
import datetime
import os
import sys

from scripts import gffnameback, gapbed, getgenepairs, gffgrep, gapseq, Rlabel, redundantfilter, findnonfunc, \
    polishtypescan, splitmodgff_genblastg, nameback, gapanno_genblastg, gapfind, modalign2origin, findID, overlapfilter, \
    gffrename


def synteny(sp1, sp2, annoType1, annoKey1, annoparentKey1, annoType2, annoKey2, annoparentKey2,
           datatype, cscore, threads):
    # Rename the input files
    script_dir, filename = os.path.split(os.path.abspath(sys.argv[0]))
    global workDir
    workDir = os.getcwd()
    sp1fa = workDir + '/' + str(sp1) + '.fa'
    sp2fa = workDir + '/' + str(sp2) + '.fa'
    sp1gff = workDir + '/' + str(sp1) + '.gff3'
    sp2gff = workDir + '/' + str(sp2) + '.gff3'
    os.chdir(workDir)

    # Perform jcvi for species1 and species2
    print('[\033[0;36mINFO\033[0m] Performing jcvi for \033[0;35m' + sp1 +
          '\033[0m and \033[0;35m' + sp2 + '\033[0m, please wait ...')
    sp1bed = workDir + '/' + str(sp1) + '.bed'
    sp2bed = workDir + '/' + str(sp2) + '.bed'
    sp1priid = workDir + '/' + str(sp1) + '.primary.id'
    sp2priid = workDir + '/' + str(sp2) + '.primary.id'
    sp1cds = workDir + '/' + str(sp1) + '.cds'
    sp2cds = workDir + '/' + str(sp2) + '.cds'
    sp1pep = workDir + '/' + str(sp1) + '.pep'
    sp2pep = workDir + '/' + str(sp2) + '.pep'
    sp1pricds = workDir + '/' + str(sp1) + '.primary.cds'
    sp2pricds = workDir + '/' + str(sp2) + '.primary.cds'
    sp1pripep = workDir + '/' + str(sp1) + '.primary.pep'
    sp2pripep = workDir + '/' + str(sp2) + '.primary.pep'
    os.system('python -m jcvi.formats.gff bed --type=' + str(annoType1) +
              ' --key=' + str(annoKey1) + ' --parent_key=' + str(annoparentKey1) +
              ' --primary_only ' + sp1gff + ' -o ' + sp1bed)
    os.system('python -m jcvi.formats.gff bed --type=' + str(annoType2) +
              ' --key=' + str(annoKey2) + ' --parent_key=' + str(annoparentKey2) +
              ' --primary_only ' + sp2gff + ' -o ' + sp2bed)
    os.system('cut -f4 ' + sp1bed + ' > ' + sp1priid)
    os.system('cut -f4 ' + sp2bed + ' > ' + sp2priid)
    os.system('gffread ' + sp1gff + ' -g ' + sp1fa + ' -x ' + sp1cds)
    os.system('seqkit grep -f ' + sp1priid + ' ' + sp1cds + ' > ' + sp1pricds)
    os.system('gffread ' + sp1gff + ' -g ' + sp1fa + ' -y ' + sp1pep + ' -S')
    os.system('seqkit grep -f ' + sp1priid + ' ' + sp1pep + ' > ' + sp1pripep)
    os.system('gffread ' + sp2gff + ' -g ' + sp2fa + ' -x ' + sp2cds)
    os.system('seqkit grep -f ' + sp2priid + ' ' + sp2cds + ' > ' + sp2pricds)
    os.system('gffread ' + sp2gff + ' -g ' + sp2fa + ' -y ' + sp2pep + ' -S')
    os.system('seqkit grep -f ' + sp2priid + ' ' + sp2pep + ' > ' + sp2pripep)
    ## perform jcvi
    global jcvi_type, jcviDir
    if str(datatype) == 'nucl':
        jcvi_type = 'jcvi_cds'
        jcviDir = workDir + '/' + jcvi_type
        os.mkdir(jcviDir)
        os.symlink(sp1pricds, jcviDir + '/' + str(sp1) + '.cds')
        os.symlink(sp1bed, jcviDir + '/' + str(sp1) + '.bed')
        os.symlink(sp2pricds, jcviDir + '/' + str(sp2) + '.cds')
        os.symlink(sp2bed, jcviDir + '/' + str(sp2) + '.bed')
    elif str(datatype) == 'prot':
        jcvi_type = 'jcvi_pep'
        jcviDir = workDir + '/' + jcvi_type
        os.mkdir(jcviDir)
        os.symlink(sp1pripep, jcviDir + '/' + str(sp1) + '.pep')
        os.symlink(sp1bed, jcviDir + '/' + str(sp1) + '.bed')
        os.symlink(sp2pripep, jcviDir + '/' + str(sp2) + '.pep')
        os.symlink(sp2bed, jcviDir + '/' + str(sp2) + '.bed')
    os.chdir(jcviDir)
    os.system('python -m jcvi.compara.catalog ortholog ' + str(sp1) + ' ' + str(sp2) +
              ' --dbtype=' + str(datatype) + ' --cscore=' + str(cscore) + ' --cpus=' +
              str(threads) + ' --no_strip_names --notex')
    print('\n[\033[0;36mINFO\033[0m]  jcvi for \033[0;35m' + sp1 + 
          '\033[0m and \033[0;35m' + sp2 + '\033[0m Done!\n')
    anchors = jcviDir + '/' + str(sp1) + '.' + str(sp2) + '.anchors'
    return anchors


if __name__ == '__main__':
    sp1fa = sys.argv[1]
    sp1gff = sys.argv[2]
    sp2fa = sys.argv[3]
    sp2gff = sys.argv[4]
    sp1 = sys.argv[5]
    sp2 = sys.argv[6]
    annoType1 = sys.argv[7]
    annoKey1 = sys.argv[8]
    annoparentKey1 = sys.argv[9]
    annoType2 = sys.argv[10]
    annoKey2 = sys.argv[11]
    annoparentKey2 = sys.argv[12]
    datatype = sys.argv[13]
    cscore = sys.argv[14]
    threads = sys.argv[15]
    synteny(sp1, sp2, annoType1, annoKey1, annoparentKey1, annoType2, annoKey2, annoparentKey2,
           datatype, cscore, threads)