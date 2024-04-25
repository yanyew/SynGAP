#!/usr/bin/env python3
# -*- coding:utf-8 -*-
import datetime
import os
import sys

from scripts import gapbed, getgenepairs, gffgrep, gapseq, gapanno_miniprot, Rlabel, redundantfilter, findnonfunc, \
    polishtypescan, splitmodgff_miniprot, gapfind, modalign2origin, findID, overlapfilter, blockR


def SynGAP(sp1, sp2, annoType1, annoKey1, annoparentKey1, annoType2, annoKey2, annoparentKey2,
           datatype, cscore, threads, kmer1, kmer2, outs, intron):
    # Rename the input files
    script_dir, filename = os.path.split(os.path.abspath(sys.argv[0]))
    global workDir
    workDir = os.getcwd()
    mission_time = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    SynGAP_Dir = workDir + '/SynGAP_' + str(sp1) + '_' + str(sp2) + '_' + mission_time
    SynGAP_workspace_Dir = SynGAP_Dir + '/workspace'
    SynGAP_results_Dir = SynGAP_Dir + '/results'
    os.mkdir(SynGAP_Dir)
    os.mkdir(SynGAP_workspace_Dir)
    os.mkdir(SynGAP_results_Dir)
    sp1fa = workDir + '/' + str(sp1) + '.fa'
    sp2fa = workDir + '/' + str(sp2) + '.fa'
    sp1gff = workDir + '/' + str(sp1) + '.gff3'
    sp2gff = workDir + '/' + str(sp2) + '.gff3'
    os.chdir(workDir)

    # Perform jcvi for species1 and species2
    print('[\033[0;36mINFO\033[0m] Performing jcvi for \033[0;35m' + sp1 +
          '\033[0m and \033[0;35m' + sp2 + '\033[0m, please wait ...')
    print('[\033[0;36mINFO\033[0m] Preparing files ...')
    ## prepare files
    os.mkdir(SynGAP_workspace_Dir + '/original')
    sp1bed = SynGAP_workspace_Dir + '/original/' + str(sp1) + '.bed'
    sp2bed = SynGAP_workspace_Dir + '/original/' + str(sp2) + '.bed'
    sp1priid = SynGAP_workspace_Dir + '/original/' + str(sp1) + '.primary.id'
    sp2priid = SynGAP_workspace_Dir + '/original/' + str(sp2) + '.primary.id'
    sp1cds = SynGAP_workspace_Dir + '/original/' + str(sp1) + '.cds'
    sp2cds = SynGAP_workspace_Dir + '/original/' + str(sp2) + '.cds'
    sp1pep = SynGAP_workspace_Dir + '/original/' + str(sp1) + '.pep'
    sp2pep = SynGAP_workspace_Dir + '/original/' + str(sp2) + '.pep'
    sp1pricds = SynGAP_workspace_Dir + '/original/' + str(sp1) + '.primary.cds'
    sp2pricds = SynGAP_workspace_Dir + '/original/' + str(sp2) + '.primary.cds'
    sp1pripep = SynGAP_workspace_Dir + '/original/' + str(sp1) + '.primary.pep'
    sp2pripep = SynGAP_workspace_Dir + '/original/' + str(sp2) + '.primary.pep'
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
        jcviDir = SynGAP_workspace_Dir + '/' + jcvi_type
        os.mkdir(jcviDir)
        os.symlink(sp1pricds, jcviDir + '/' + str(sp1) + '.cds')
        os.symlink(sp1bed, jcviDir + '/' + str(sp1) + '.bed')
        os.symlink(sp2pricds, jcviDir + '/' + str(sp2) + '.cds')
        os.symlink(sp2bed, jcviDir + '/' + str(sp2) + '.bed')
    elif str(datatype) == 'prot':
        jcvi_type = 'jcvi_pep'
        jcviDir = SynGAP_workspace_Dir + '/' + jcvi_type
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

    # Find gaps from the .anchors file
    os.chdir(SynGAP_workspace_Dir)
    print('[\033[0;36mINFO\033[0m] Finding gaps from the .anchors file, please wait ...')
    anchors = SynGAP_workspace_Dir + '/' + str(sp1) + '.' + str(sp2) + '.anchors'
    anchors_gap = SynGAP_workspace_Dir + '/' + str(sp1) + '.' + str(sp2) + '.anchors.gap'
    os.symlink(jcviDir + '/' + str(sp1) + '.' + str(sp2) + '.anchors', anchors)
    print('[\033[0;36mINFO\033[0m] Finding ...')
    gapfind.gap_find(sp1bed, sp2bed, anchors, anchors_gap)
    print('[\033[0;36mINFO\033[0m] Running Done!\n')
    print('[\033[0;36mINFO\033[0m] Extracting syntenic gene pairs from jcvi results for '
          '\033[0;35m' + str(sp1) + '\033[0m and \033[0;35m' + str(sp2) + '\033[0m , please wait ...')
    Sgenepairs = SynGAP_workspace_Dir + '/' + str(sp1) + '.' + str(sp2) + '.anchors.genepairs'
    getgenepairs.get_genepair(anchors, Sgenepairs)
    print('[\033[0;36mINFO\033[0m] Running Done!\n')
    os.symlink(anchors, SynGAP_results_Dir + '/' + str(sp1) + '.' + str(sp2) + '.anchors')
    os.symlink(anchors_gap, SynGAP_results_Dir + '/' + str(sp1) + '.' + str(sp2) + '.anchors.gap')
    os.symlink(Sgenepairs, SynGAP_results_Dir + '/' + str(sp1) + '.' + str(sp2) + '.anchors.genepairs')

    # Calculate R value for blocks in .anchors
    print('[\033[0;36mINFO\033[0m] Calculate R value for blocks in .anchors, please wait ...')
    outR = SynGAP_workspace_Dir + '/' + str(sp1) + '.' + str(sp2) + '.anchors.R.xls'
    sp1blockRbed = SynGAP_workspace_Dir + '/' + str(sp1) + '.block.bed'
    sp2blockRbed = SynGAP_workspace_Dir + '/' + str(sp2) + '.block.bed'
    blockR.blockR(anchors, sp1blockRbed, sp2blockRbed, sp1pripep, sp2pripep, str(threads), sp1bed, sp2bed, outR)

    # Extract .bed information for the gaps from .anchors.gap
    print('[\033[0;36mINFO\033[0m] Extracting .bed information for the gaps, please wait ...')
    sp1gapbed = SynGAP_workspace_Dir + '/' + str(sp1) + '.gap.bed'
    sp2gapbed = SynGAP_workspace_Dir + '/' + str(sp2) + '.gap.bed'
    gapbed.gapbed(anchors_gap, sp1bed, sp2bed, sp1gapbed, sp2gapbed)
    print('[\033[0;36mINFO\033[0m] Running Done!\n')

    # Prepare sequences data for gap annotation
    print('[\033[0;36mINFO\033[0m] Preparing sequences for gap annotation, please wait ...')
    os.mkdir(SynGAP_workspace_Dir + '/gapanno')
    os.chdir(SynGAP_workspace_Dir + '/gapanno')
    gapseq.gapseq(sp1fa, sp2fa, sp1pep, sp2pep, sp1gapbed, sp2gapbed)
    print('[\033[0;36mINFO\033[0m] Running Done!\n')

    # Perform gap annotation
    print('[\033[0;36mINFO\033[0m] Performing gap annotation, please wait ...')
    gapanno_miniprot.gapanno(str(script_dir) + '/bin/miniprot', str(kmer1), str(kmer2), str(outs), str(intron),
                             str(threads))
    print('[\033[0;36mINFO\033[0m] Running Done!\n')

    # Split the modified annotation file
    print('[\033[0;36mINFO\033[0m] Spliting the raw modified annotations file, please wait ...')
    os.chdir(SynGAP_workspace_Dir)
    modgff = SynGAP_workspace_Dir + '/' + 'modified.raw.gff'
    sp1modgff = SynGAP_workspace_Dir + '/' + str(sp1) + '.modified.gff'
    sp2modgff = SynGAP_workspace_Dir + '/' + str(sp2) + '.modified.gff'
    splitmodgff_miniprot.split(str(sp1), str(sp2), sp1gapbed, sp2gapbed, modgff, sp1modgff, sp2modgff)
    print('[\033[0;36mINFO\033[0m] Running Done!\n')

    # Filter out the modified annotations that are redundant and overlapping with the annotations from genome .gff3
    sp1modbed = SynGAP_workspace_Dir + '/filter/' + str(sp1) + '.modified.bed'
    sp2modbed = SynGAP_workspace_Dir + '/filter/' + str(sp2) + '.modified.bed'
    sp1modgff_disredundant = SynGAP_workspace_Dir + '/filter/' + str(sp1) + '.modified.disredundant.gff'
    sp2modgff_disredundant = SynGAP_workspace_Dir + '/filter/' + str(sp2) + '.modified.disredundant.gff'
    sp1modgff_filtered = SynGAP_workspace_Dir + '/' + str(sp1) + '.modified.filtered.gff'
    sp2modgff_filtered = SynGAP_workspace_Dir + '/' + str(sp2) + '.modified.filtered.gff'
    print('[\033[0;36mINFO\033[0m] Filtering out the redundant modified annotations, please wait ...')
    os.mkdir(SynGAP_workspace_Dir + '/filter')
    os.system(
        'python -m jcvi.formats.gff bed --type=mRNA --key=ID --parent_key=Parent ' + sp1modgff + ' -o ' + sp1modbed)
    os.system(
        'python -m jcvi.formats.gff bed --type=mRNA --key=ID --parent_key=Parent ' + sp2modgff + ' -o ' + sp2modbed)
    redundantfilter.redundantfilter(sp1modgff, sp1modbed, sp1modgff_disredundant)
    print('[\033[0;36mINFO\033[0m] Running Done!')
    redundantfilter.redundantfilter(sp2modgff, sp2modbed, sp2modgff_disredundant)
    print('[\033[0;36mINFO\033[0m] Running Done!')
    print('[\033[0;36mINFO\033[0m] Filtering out the modified annotations overlapping with the '
          'original genome annotations, please wait ...')
    overlapfilter.overlapfilter(
        sp1gff, str(annoType1), str(annoKey1), sp1modgff_disredundant, sp1modgff_filtered)
    print('[\033[0;36mINFO\033[0m] Running Done!')
    overlapfilter.overlapfilter(
        sp2gff, str(annoType2), str(annoKey2), sp2modgff_disredundant, sp2modgff_filtered)
    print('[\033[0;36mINFO\033[0m] Running Done!\n')

    # Calculate R value for modified annotations
    print('[\033[0;36mINFO\033[0m] Labeling \033[0;33mR value\033[0m for modified annotations, please wait ...')
    os.mkdir(SynGAP_workspace_Dir + '/Rlabel')
    os.chdir(SynGAP_workspace_Dir + '/Rlabel')
    sp1modpep = SynGAP_workspace_Dir + '/Rlabel/' + str(sp1) + '.filtered.pep'
    sp2modpep = SynGAP_workspace_Dir + '/Rlabel/' + str(sp2) + '.filtered.pep'
    os.system('gffread ' + sp1modgff_filtered + ' -g ' + sp1fa + ' -y ' + sp1modpep + ' -S')
    os.system('gffread ' + sp2modgff_filtered + ' -g ' + sp2fa + ' -y ' + sp2modpep + ' -S')
    sp1modpep_originalid = SynGAP_workspace_Dir + '/Rlabel/' + str(sp1) + '.filtered.origin.id'
    sp2modpep_originalid = SynGAP_workspace_Dir + '/Rlabel/' + str(sp2) + '.filtered.origin.id'
    sp1modpep_originpep = SynGAP_workspace_Dir + '/Rlabel/' + str(sp1) + '.filtered.origin.pep'
    sp2modpep_originpep = SynGAP_workspace_Dir + '/Rlabel/' + str(sp2) + '.filtered.origin.pep'
    os.system("grep mRNA " + sp1modgff_filtered + " | cut -f9 | cut -d';' -f1 | sed 's/ID=" + str(sp1) +
              "_//g' | sed 's/_[1-9][0-9]$//g' | sed 's/_[1-9]$//g' | sort | uniq > " + sp1modpep_originalid)
    os.system("grep mRNA " + sp2modgff_filtered + " | cut -f9 | cut -d';' -f1 | sed 's/ID=" + str(sp2) +
              "_//g' | sed 's/_[1-9][0-9]$//g' | sed 's/_[1-9]$//g' | sort | uniq > " + sp2modpep_originalid)
    os.system("seqkit grep -f " + sp1modpep_originalid + " " + sp2pripep + " > " + sp1modpep_originpep)
    os.system("seqkit grep -f " + sp2modpep_originalid + " " + sp1pripep + " > " + sp2modpep_originpep)
    sp1R = SynGAP_workspace_Dir + '/' + str(sp1) + '.modified.filtered.R.xls'
    sp2R = SynGAP_workspace_Dir + '/' + str(sp2) + '.modified.filtered.R.xls'
    sp1modgff_R = SynGAP_workspace_Dir + '/' + str(sp1) + '.modified.filtered.R.gff'
    sp2modgff_R = SynGAP_workspace_Dir + '/' + str(sp2) + '.modified.filtered.R.gff'
    print('[\033[0;36mINFO\033[0m] Performing global alignment of '
          'modified annotations and original annotations for ' + str(sp1) + ', please wait ...')
    print('[\033[0;36mINFO\033[0m] Load file `' + sp1modpep + '`')
    print('[\033[0;36mINFO\033[0m] Load file `' + sp1modpep_originpep + '`')
    modalign2origin.needle(str(sp1), sp1modpep, sp1modpep_originpep, str(threads))
    print('[\033[0;36mINFO\033[0m] Running Done!')
    print('[\033[0;36mINFO\033[0m] Performing global alignment of '
          'modified annotations and original annotations for ' + str(sp2) + ', please wait ...')
    print('[\033[0;36mINFO\033[0m] Load file `' + sp2modpep + '`')
    print('[\033[0;36mINFO\033[0m] Load file `' + sp2modpep_originpep + '`')
    modalign2origin.needle(str(sp2), sp2modpep, sp2modpep_originpep, str(threads))
    print('[\033[0;36mINFO\033[0m] Running Done!')
    print('[\033[0;36mINFO\033[0m] Calculating and labeling \033[0;33mR value\033[0m '
          'for modified annotations of \033[0;35m' + str(sp1) + '\033[0m, please wait ...')
    Rlabel.Rlabel(str(sp1), sp1blockRbed, sp2bed, sp2blockRbed, sp1modgff_filtered, sp1R, sp1modgff_R)
    print('[\033[0;36mINFO\033[0m] Running Done!')
    print('[\033[0;36mINFO\033[0m] Calculating and labeling \033[0;33mR value\033[0m '
          'for modified annotations of \033[0;35m' + str(sp2) + '\033[0m, please wait ...')
    Rlabel.Rlabel(str(sp2), sp2blockRbed, sp1bed, sp1blockRbed, sp2modgff_filtered, sp2R, sp2modgff_R)
    print('[\033[0;36mINFO\033[0m] Running Done!\n')

    # Scan out the modified annotations that are not overlapped with the annotations from genome .gff3
    sp1modgff_miss_annotated = SynGAP_workspace_Dir + '/' + str(sp1) + '.modified.filtered.R.miss_annotated.gff'
    sp1modgff_mis_annotated = SynGAP_workspace_Dir + '/' + str(sp1) + '.modified.filtered.R.mis_annotated.gff'
    sp2modgff_miss_annotated = SynGAP_workspace_Dir + '/' + str(sp2) + '.modified.filtered.R.miss_annotated.gff'
    sp2modgff_mis_annotated = SynGAP_workspace_Dir + '/' + str(sp2) + '.modified.filtered.R.mis_annotated.gff'
    print('[\033[0;36mINFO\033[0m] Scanning the modified annotations not overlapping with the '
          'original genome annotations, please wait ...')
    polishtypescan.polishtypescan(sp1modgff_R, sp1modgff_miss_annotated, sp1modgff_mis_annotated)
    print('[\033[0;36mINFO\033[0m] Running Done!')
    polishtypescan.polishtypescan(sp2modgff_R, sp2modgff_miss_annotated, sp2modgff_mis_annotated)
    print('[\033[0;36mINFO\033[0m] Running Done!')

    # Perform the CPC to predict coding posibility for original genome annotations
    print('[\033[0;36mINFO\033[0m] Performing CPC to predict coding posibility for '
          'original genome annotations, please wait ...')
    os.mkdir(SynGAP_workspace_Dir + '/CPC')
    sp1cpc = SynGAP_workspace_Dir + '/CPC/' + str(sp1) + '.cpc2output'
    sp2cpc = SynGAP_workspace_Dir + '/CPC/' + str(sp2) + '.cpc2output'
    os.system('python ' + str(script_dir) + '/bin/CPC2_standalone-1.0.1/bin/CPC2.py -i ' + sp1pricds + ' -o ' + sp1cpc)
    os.system('python ' + str(script_dir) + '/bin/CPC2_standalone-1.0.1/bin/CPC2.py -i ' + sp2pricds + ' -o ' + sp2cpc)
    sp1id_coding = SynGAP_workspace_Dir + '/CPC/' + str(sp1) + '.coding.id'
    sp1id_noncoding = SynGAP_workspace_Dir + '/CPC/' + str(sp1) + '.noncoding.id'
    sp2id_coding = SynGAP_workspace_Dir + '/CPC/' + str(sp2) + '.coding.id'
    sp2id_noncoding = SynGAP_workspace_Dir + '/CPC/' + str(sp2) + '.noncoding.id'
    os.system("awk '{if($7<=0.5) print $1}' " + sp1cpc + ".txt | cut -f1 | sort | uniq > " + sp1id_noncoding)
    os.system("awk '{if($7>0.5) print $1}' " + sp1cpc + ".txt | cut -f1 | sort | uniq > " + sp1id_coding)
    os.system("awk '{if($7<=0.5) print $1}' " + sp2cpc + ".txt | cut -f1 | sort | uniq > " + sp2id_noncoding)
    os.system("awk '{if($7>0.5) print $1}' " + sp2cpc + ".txt | cut -f1 | sort | uniq > " + sp2id_coding)
    print('[\033[0;36mINFO\033[0m] Running Done!\n')

    # Perform gene function prediction for original genome annotations
    print(
        '[\033[0;36mINFO\033[0m] Performing gene function prediction for original genome annotations, please wait ...')
    os.mkdir(SynGAP_workspace_Dir + '/function')
    sp1func = SynGAP_workspace_Dir + '/function/' + str(sp1) + '2uniprot.txt'
    sp2func = SynGAP_workspace_Dir + '/function/' + str(sp2) + '2uniprot.txt'
    sp1id_function = SynGAP_workspace_Dir + '/function/' + str(sp1) + '.function.id'
    sp1id_nonfunction = SynGAP_workspace_Dir + '/function/' + str(sp1) + '.nonfunction.id'
    sp2id_function = SynGAP_workspace_Dir + '/function/' + str(sp2) + '.function.id'
    sp2id_nonfunctionn = SynGAP_workspace_Dir + '/function/' + str(sp2) + '.nonfunction.id'
    os.system('diamond blastp --db ' + str(script_dir) + '/bin/SwissProt/uniprot_sprot.fasta' + ' --query ' +
              sp1pripep + ' -e 1e-8 --outfmt 6 --max-target-seqs 5 --threads ' + str(threads) +
              ' --quiet --out ' + sp1func)
    os.system('diamond blastp --db ' + str(script_dir) + '/bin/SwissProt/uniprot_sprot.fasta' + ' --query ' +
              sp2pripep + ' -e 1e-8 --outfmt 6 --max-target-seqs 5 --threads ' + str(threads) +
              ' --quiet --out ' + sp2func)
    os.system('cut -f1 ' + sp1func + ' | sort | uniq > ' + sp1id_function)
    os.system('cut -f1 ' + sp2func + ' | sort | uniq > ' + sp2id_function)
    findnonfunc.filter(sp1id_function, sp1priid, sp1id_nonfunction)
    findnonfunc.filter(sp2id_function, sp2priid, sp2id_nonfunctionn)
    print('[\033[0;36mINFO\033[0m] Running Done!\n')

    # Classify the modified annotations based on coding posibility and predicted gene function
    print('[\033[0;36mINFO\033[0m] Classifying the modified annotations based on coding posibility '
          'and predicted gene function of corresponding original annotations, please wait ...')
    os.mkdir(SynGAP_workspace_Dir + '/cluster')
    sp1id_I = SynGAP_workspace_Dir + '/cluster/' + str(sp1) + '.I.id'
    sp1id_IIa = SynGAP_workspace_Dir + '/cluster/' + str(sp1) + '.IIa.id'
    sp1id_IIb = SynGAP_workspace_Dir + '/cluster/' + str(sp1) + '.IIb.id'
    sp1id_III = SynGAP_workspace_Dir + '/cluster/' + str(sp1) + '.III.id'
    sp2id_I = SynGAP_workspace_Dir + '/cluster/' + str(sp2) + '.I.id'
    sp2id_IIa = SynGAP_workspace_Dir + '/cluster/' + str(sp2) + '.IIa.id'
    sp2id_IIb = SynGAP_workspace_Dir + '/cluster/' + str(sp2) + '.IIb.id'
    sp2id_III = SynGAP_workspace_Dir + '/cluster/' + str(sp2) + '.III.id'
    findID.filter(sp1id_coding, sp1id_function, sp1id_I)
    findID.filter(sp1id_noncoding, sp1id_nonfunction, sp1id_III)
    findID.filter(sp1id_coding, sp1id_nonfunction, sp1id_IIa)
    findID.filter(sp1id_noncoding, sp1id_function, sp1id_IIb)
    findID.filter(sp2id_coding, sp2id_function, sp2id_I)
    findID.filter(sp2id_noncoding, sp2id_nonfunctionn, sp2id_III)
    findID.filter(sp2id_coding, sp2id_nonfunctionn, sp2id_IIa)
    findID.filter(sp2id_noncoding, sp2id_function, sp2id_IIb)
    sp1modgff_I = SynGAP_workspace_Dir + '/cluster/' + str(sp1) + '.modified.filtered.R.I.gff'
    sp1modgff_III = SynGAP_workspace_Dir + '/cluster/' + str(sp1) + '.modified.filtered.R.III.gff'
    sp1modgff_IIa = SynGAP_workspace_Dir + '/cluster/' + str(sp1) + '.modified.filtered.R.IIa.gff'
    sp1modgff_IIb = SynGAP_workspace_Dir + '/cluster/' + str(sp1) + '.modified.filtered.R.IIb.gff'
    sp2modgff_I = SynGAP_workspace_Dir + '/cluster/' + str(sp2) + '.modified.filtered.R.I.gff'
    sp2modgff_III = SynGAP_workspace_Dir + '/cluster/' + str(sp2) + '.modified.filtered.R.III.gff'
    sp2modgff_IIa = SynGAP_workspace_Dir + '/cluster/' + str(sp2) + '.modified.filtered.R.IIa.gff'
    sp2modgff_IIb = SynGAP_workspace_Dir + '/cluster/' + str(sp2) + '.modified.filtered.R.IIb.gff'
    gffgrep.gffgrep(sp1modgff_R, sp2id_I, sp2id_IIa, sp2id_IIb, sp2id_III,
                    sp1modgff_I, sp1modgff_IIa, sp1modgff_IIb, sp1modgff_III)
    gffgrep.gffgrep(sp2modgff_R, sp1id_I, sp1id_IIa, sp1id_IIb, sp1id_III,
                    sp2modgff_I, sp2modgff_IIa, sp2modgff_IIb, sp2modgff_III)
    sp1modgff_miss_annotated_I = SynGAP_workspace_Dir + '/cluster/' + str(
        sp1) + '.modified.filtered.R.miss_annotated.I.gff'
    sp1modgff_mis_annotated_I = SynGAP_workspace_Dir + '/cluster/' + str(
        sp1) + '.modified.filtered.R.mis_annotated.I.gff'
    sp1modgff_miss_annotated_III = SynGAP_workspace_Dir + '/cluster/' + str(
        sp1) + '.modified.filtered.R.miss_annotated.III.gff'
    sp1modgff_mis_annotated_III = SynGAP_workspace_Dir + '/cluster/' + str(
        sp1) + '.modified.filtered.R.mis_annotated.III.gff'
    sp1modgff_miss_annotated_IIa = SynGAP_workspace_Dir + '/cluster/' + str(
        sp1) + '.modified.filtered.R.miss_annotated.IIa.gff'
    sp1modgff_mis_annotated_IIa = SynGAP_workspace_Dir + '/cluster/' + str(
        sp1) + '.modified.filtered.R.mis_annotated.IIa.gff'
    sp1modgff_miss_annotated_IIb = SynGAP_workspace_Dir + '/cluster/' + str(
        sp1) + '.modified.filtered.R.miss_annotated.IIb.gff'
    sp1modgff_mis_annotated_IIb = SynGAP_workspace_Dir + '/cluster/' + str(
        sp1) + '.modified.filtered.R.mis_annotated.IIb.gff'
    sp2modgff_miss_annotated_I = SynGAP_workspace_Dir + '/cluster/' + str(
        sp2) + '.modified.filtered.R.miss_annotated.I.gff'
    sp2modgff_mis_annotated_I = SynGAP_workspace_Dir + '/cluster/' + str(
        sp2) + '.modified.filtered.R.mis_annotated.I.gff'
    sp2modgff_miss_annotated_III = SynGAP_workspace_Dir + '/cluster/' + str(
        sp2) + '.modified.filtered.R.miss_annotated.III.gff'
    sp2modgff_mis_annotated_III = SynGAP_workspace_Dir + '/cluster/' + str(
        sp2) + '.modified.filtered.R.mis_annotated.III.gff'
    sp2modgff_miss_annotated_IIa = SynGAP_workspace_Dir + '/cluster/' + str(
        sp2) + '.modified.filtered.R.miss_annotated.IIa.gff'
    sp2modgff_mis_annotated_IIa = SynGAP_workspace_Dir + '/cluster/' + str(
        sp2) + '.modified.filtered.R.mis_annotated.IIa.gff'
    sp2modgff_miss_annotated_IIb = SynGAP_workspace_Dir + '/cluster/' + str(
        sp2) + '.modified.filtered.R.miss_annotated.IIb.gff'
    sp2modgff_mis_annotated_IIb = SynGAP_workspace_Dir + '/cluster/' + str(
        sp2) + '.modified.filtered.R.mis_annotated.IIb.gff'
    polishtypescan.polishtypescan(sp1modgff_I, sp1modgff_miss_annotated_I, sp1modgff_mis_annotated_I)
    polishtypescan.polishtypescan(sp1modgff_III, sp1modgff_miss_annotated_III, sp1modgff_mis_annotated_III)
    polishtypescan.polishtypescan(sp1modgff_IIa, sp1modgff_miss_annotated_IIa, sp1modgff_mis_annotated_IIa)
    polishtypescan.polishtypescan(sp1modgff_IIb, sp1modgff_miss_annotated_IIb, sp1modgff_mis_annotated_IIb)
    polishtypescan.polishtypescan(sp2modgff_I, sp2modgff_miss_annotated_I, sp2modgff_mis_annotated_I)
    polishtypescan.polishtypescan(sp2modgff_III, sp2modgff_miss_annotated_III, sp2modgff_mis_annotated_III)
    polishtypescan.polishtypescan(sp2modgff_IIa, sp2modgff_miss_annotated_IIa, sp2modgff_mis_annotated_IIa)
    polishtypescan.polishtypescan(sp2modgff_IIb, sp2modgff_miss_annotated_IIb, sp2modgff_mis_annotated_IIb)
    sp1GAPgff = SynGAP_results_Dir + '/' + str(sp1) + '.SynGAP.gff3'
    sp2GAPgff = SynGAP_results_Dir + '/' + str(sp2) + '.SynGAP.gff3'
    sp1GAPcleangff = SynGAP_results_Dir + '/' + str(sp1) + '.SynGAP.clean.gff3'
    sp2GAPcleangff = SynGAP_results_Dir + '/' + str(sp2) + '.SynGAP.clean.gff3'
    sp1GAPcleanmiss_annotatedgff = SynGAP_results_Dir + '/' + str(sp1) + '.SynGAP.clean.miss_annotated.gff3'
    sp1GAPcleanmis_annotatedgff = SynGAP_results_Dir + '/' + str(sp1) + '.SynGAP.clean.mis_annotated.gff3'
    sp2GAPcleanmiss_annotatedgff = SynGAP_results_Dir + '/' + str(sp2) + '.SynGAP.clean.miss_annotated.gff3'
    sp2GAPcleanmis_annotatedgff = SynGAP_results_Dir + '/' + str(sp2) + '.SynGAP.clean.mis_annotated.gff3'
    os.chdir(SynGAP_results_Dir)
    os.symlink(sp1modgff_I, sp1GAPcleangff)
    os.symlink(sp2modgff_I, sp2GAPcleangff)
    os.symlink(sp1modgff_miss_annotated_I, sp1GAPcleanmiss_annotatedgff)
    os.symlink(sp1modgff_mis_annotated_I, sp1GAPcleanmis_annotatedgff)
    os.symlink(sp2modgff_miss_annotated_I, sp2GAPcleanmiss_annotatedgff)
    os.symlink(sp2modgff_mis_annotated_I, sp2GAPcleanmis_annotatedgff)
    os.system('cat ' + sp1gff + ' ' + sp1modgff_I + ' > ' + sp1GAPgff)
    os.system('cat ' + sp2gff + ' ' + sp2modgff_I + ' > ' + sp2GAPgff)
    print('[\033[0;36mINFO\033[0m] Running Done!\n')

    # Perform jcvi for species1 and species2 after SynGAP
    print('[\033[0;36mINFO\033[0m] Performing jcvi for \033[0;35m' + str(sp1) +
          '\033[0m and \033[0;35m' + str(sp2) + '\033[0m after SynGAP, please wait ...')
    print('[\033[0;36mINFO\033[0m] Preparing files ...')
    jcviafterSynGAP_Dir = SynGAP_workspace_Dir + '/jcviafterSynGAP'
    os.mkdir(jcviafterSynGAP_Dir)
    sp1GAPbed = jcviafterSynGAP_Dir + '/' + str(sp1) + '.SynGAP.bed'
    sp2GAPbed = jcviafterSynGAP_Dir + '/' + str(sp2) + '.SynGAP.bed'
    sp1GAPpriid = jcviafterSynGAP_Dir + '/' + str(sp1) + '.SynGAP.primary.id'
    sp2GAPpriid = jcviafterSynGAP_Dir + '/' + str(sp2) + '.SynGAP.primary.id'
    sp1GAPcds = jcviafterSynGAP_Dir + '/' + str(sp1) + '.SynGAP.cds'
    sp2GAPcds = jcviafterSynGAP_Dir + '/' + str(sp2) + '.SynGAP.cds'
    sp1GAPpep = jcviafterSynGAP_Dir + '/' + str(sp1) + '.SynGAP.pep'
    sp2GAPpep = jcviafterSynGAP_Dir + '/' + str(sp2) + '.SynGAP.pep'
    sp1GAPpricds = jcviafterSynGAP_Dir + '/' + str(sp1) + '.SynGAP.primary.cds'
    sp2GAPpricds = jcviafterSynGAP_Dir + '/' + str(sp2) + '.SynGAP.primary.cds'
    sp1GAPpripep = jcviafterSynGAP_Dir + '/' + str(sp1) + '.SynGAP.primary.pep'
    sp2GAPpripep = jcviafterSynGAP_Dir + '/' + str(sp2) + '.SynGAP.primary.pep'
    os.system('python -m jcvi.formats.gff bed --type=' + str(annoType1) +
              ' --key=' + str(annoKey1) + ' --parent_key=' + str(annoparentKey1) +
              ' --primary_only ' + sp1GAPgff + ' -o ' + sp1GAPbed)
    os.system('python -m jcvi.formats.gff bed --type=' + str(annoType2) +
              ' --key=' + str(annoKey2) + ' --parent_key=' + str(annoparentKey2) +
              ' --primary_only ' + sp2GAPgff + ' -o ' + sp2GAPbed)
    os.system('cut -f4 ' + sp1GAPbed + ' > ' + sp1GAPpriid)
    os.system('cut -f4 ' + sp2GAPbed + ' > ' + sp2GAPpriid)
    os.system('gffread ' + sp1GAPgff + ' -g ' + sp1fa + ' -x ' + sp1GAPcds)
    os.system('seqkit grep -f ' + sp1GAPpriid + ' ' + sp1GAPcds + ' > ' + sp1GAPpricds)
    os.system('gffread ' + sp1GAPgff + ' -g ' + sp1fa + ' -y ' + sp1GAPpep + ' -S')
    os.system('seqkit grep -f ' + sp1GAPpriid + ' ' + sp1GAPpep + ' > ' + sp1GAPpripep)
    os.system('gffread ' + sp2GAPgff + ' -g ' + sp2fa + ' -x ' + sp2GAPcds)
    os.system('seqkit grep -f ' + sp2GAPpriid + ' ' + sp2GAPcds + ' > ' + sp2GAPpricds)
    os.system('gffread ' + sp2GAPgff + ' -g ' + sp2fa + ' -y ' + sp2GAPpep + ' -S')
    os.system('seqkit grep -f ' + sp2GAPpriid + ' ' + sp2GAPpep + ' > ' + sp2GAPpripep)
    global jcviafterSynGAPDir
    if str(datatype) == 'nucl':
        jcvi_type = 'jcvi_cds'
        jcviafterSynGAPDir = jcviafterSynGAP_Dir + '/' + jcvi_type
        os.mkdir(jcviafterSynGAPDir)
        os.symlink(sp1GAPpricds, jcviafterSynGAPDir + '/' + str(sp1) + '.cds')
        os.symlink(sp1GAPbed, jcviafterSynGAPDir + '/' + str(sp1) + '.bed')
        os.symlink(sp2GAPpricds, jcviafterSynGAPDir + '/' + str(sp2) + '.cds')
        os.symlink(sp2GAPbed, jcviafterSynGAPDir + '/' + str(sp2) + '.bed')
    elif str(datatype) == 'prot':
        jcvi_type = 'jcvi_pep'
        jcviafterSynGAPDir = jcviafterSynGAP_Dir + '/' + jcvi_type
        os.mkdir(jcviafterSynGAPDir)
        os.symlink(sp1GAPpripep, jcviafterSynGAPDir + '/' + str(sp1) + '.pep')
        os.symlink(sp1GAPbed, jcviafterSynGAPDir + '/' + str(sp1) + '.bed')
        os.symlink(sp2GAPpripep, jcviafterSynGAPDir + '/' + str(sp2) + '.pep')
        os.symlink(sp2GAPbed, jcviafterSynGAPDir + '/' + str(sp2) + '.bed')
    os.chdir(jcviafterSynGAPDir)
    os.system('python -m jcvi.compara.catalog ortholog ' + str(sp1) + ' ' + str(sp2) +
              ' --dbtype=' + str(datatype) + ' --cscore=' + str(cscore) + ' --cpus=' +
              str(threads) + ' --no_strip_names --notex')
    print('\n[\033[0;36mINFO\033[0m] jcvi for \033[0;35m' + str(sp1) +
          '\033[0m and \033[0;35m' + str(sp2) + '\033[0m after SynGAP Done!\n')

    # Extract syntenic gene pairs from jcvi results for species1 and species2 after SynGAP
    anchors_GAP = jcviafterSynGAP_Dir + '/' + jcvi_type + '/' + str(
        sp1) + '.' + str(sp2) + '.anchors'
    Sgenepairs_GAP = jcviafterSynGAP_Dir + '/' + jcvi_type + '/' + str(
        sp1) + '.' + str(sp2) + '.anchors.SynGAP.genepairs'
    print('[\033[0;36mINFO\033[0m] Extracting syntenic gene pairs from jcvi results for '
          '\033[0;35m' + str(sp1) + '\033[0m and \033[0;35m' + str(sp2) + '\033[0m after SynGAP, please wait ...')
    getgenepairs.get_genepair(anchors_GAP, Sgenepairs_GAP)
    os.chdir(SynGAP_results_Dir)
    os.symlink(Sgenepairs_GAP, str(sp1) + '.' + str(sp2) + '.anchors.SynGAP.genepairs')
    print('[\033[0;36mINFO\033[0m] Running Done!\n')

    # Finished message
    print('[\033[0;36mINFO\033[0m] SynGAP analysis for \033[0;35m' + str(sp1) +
          '\033[0m and \033[0;35m' + str(sp2) + '\033[0m Done!')
    print('[\033[0;36mINFO\033[0m] Please check the result files in `\033[0;35m' + SynGAP_results_Dir + '\033[0m`\n')
    os.chdir(workDir)
    return SynGAP_results_Dir


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
    kmer1 = sys.argv[16]
    kmer2 = sys.argv[17]
    outs = sys.argv[18]
    intron = sys.argv[19]
    SynGAP(sp1, sp2, annoType1, annoKey1, annoparentKey1, annoType2, annoKey2, annoparentKey2,
           datatype, cscore, threads, kmer1, kmer2, outs, intron)