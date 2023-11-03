#!/usr/bin/env python3
# -*- coding:utf-8 -*-
import os
import sys
import datetime
from scripts import gffnameback, gapbed, getgenepairs, gffgrep, gapseq, Rlabel, redundantfilter, findnonfunc, \
    polishtypescan, splitmodgff_genblastg, nameback, gapanno_genblastg, gapfind, modalign2origin, findID, overlapfilter, \
    gffrename


def SynGAP(sp1, sp2, annoType1, annoKey1, annoparentKey1, annoType2, annoKey2, annoparentKey2,
           datatype, cscore, threads, evalue, rank, coverage):
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

    # Rename the gene_ids in annotation files
    sp1gff_rename = SynGAP_workspace_Dir + '/renamed_id/' + str(sp1) + '.renamed.gff3'
    sp2gff_rename = SynGAP_workspace_Dir + '/renamed_id/' + str(sp2) + '.renamed.gff3'
    sp1mRNAmap = SynGAP_workspace_Dir + '/renamed_id/' + str(sp1) + '.rename.mRNAmap'
    sp2mRNAmap = SynGAP_workspace_Dir + '/renamed_id/' + str(sp2) + '.rename.mRNAmap'
    os.mkdir(SynGAP_workspace_Dir + '/original_id')
    os.mkdir(SynGAP_workspace_Dir + '/renamed_id')
    print('[\033[0;36mINFO\033[0m] Renaming the gene ids in annotation files, please wait ...')
    gffrename.rename_gff(sp1gff, str(sp1), str(annoType1), str(annoKey1),
                         str(annoparentKey1), sp1gff_rename, sp1mRNAmap)
    print('[\033[0;36mINFO\033[0m] Running Done!')
    gffrename.rename_gff(sp2gff, str(sp2), str(annoType2), str(annoKey2),
                         str(annoparentKey2), sp2gff_rename, sp2mRNAmap)
    print('[\033[0;36mINFO\033[0m] Running Done!\n')

    # Perform MCScan for species1 and species2
    print('[\033[0;36mINFO\033[0m] Performing MCScan for \033[0;35m' + sp1 + '\033[0m '
          'and \033[0;35m' + sp2 + '\033[0m, please wait ...')
    print('[\033[0;36mINFO\033[0m] Preparing files ...')
    ## prepare files for original_id
    sp1bed = SynGAP_workspace_Dir + '/original_id/' + str(sp1) + '.bed'
    sp2bed = SynGAP_workspace_Dir + '/original_id/' + str(sp2) + '.bed'
    sp1priid = SynGAP_workspace_Dir + '/original_id/' + str(sp1) + '.primary.id'
    sp2priid = SynGAP_workspace_Dir + '/original_id/' + str(sp2) + '.primary.id'
    sp1cds = SynGAP_workspace_Dir + '/original_id/' + str(sp1) + '.cds'
    sp2cds = SynGAP_workspace_Dir + '/original_id/' + str(sp2) + '.cds'
    sp1pep = SynGAP_workspace_Dir + '/original_id/' + str(sp1) + '.pep'
    sp2pep = SynGAP_workspace_Dir + '/original_id/' + str(sp2) + '.pep'
    sp1pricds = SynGAP_workspace_Dir + '/original_id/' + str(sp1) + '.primary.cds'
    sp2pricds = SynGAP_workspace_Dir + '/original_id/' + str(sp2) + '.primary.cds'
    sp1pripep = SynGAP_workspace_Dir + '/original_id/' + str(sp1) + '.primary.pep'
    sp2pripep = SynGAP_workspace_Dir + '/original_id/' + str(sp2) + '.primary.pep'
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
    ## prepare files for renamed_id
    sp1bed_rename = SynGAP_workspace_Dir + '/renamed_id/' + str(sp1) + '.rename.bed'
    sp2bed_rename = SynGAP_workspace_Dir + '/renamed_id/' + str(sp2) + '.rename.bed'
    sp1priid_rename = SynGAP_workspace_Dir + '/renamed_id/' + str(sp1) + '.rename.primary.id'
    sp2priid_rename = SynGAP_workspace_Dir + '/renamed_id/' + str(sp2) + '.rename.primary.id'
    sp1cds_rename = SynGAP_workspace_Dir + '/renamed_id/' + str(sp1) + '.rename.cds'
    sp2cds_rename = SynGAP_workspace_Dir + '/renamed_id/' + str(sp2) + '.rename.cds'
    sp1pep_rename = SynGAP_workspace_Dir + '/renamed_id/' + str(sp1) + '.rename.pep'
    sp2pep_rename = SynGAP_workspace_Dir + '/renamed_id/' + str(sp2) + '.rename.pep'
    sp1pricds_rename = SynGAP_workspace_Dir + '/renamed_id/' + str(sp1) + '.rename.primary.cds'
    sp2pricds_rename = SynGAP_workspace_Dir + '/renamed_id/' + str(sp2) + '.rename.primary.cds'
    sp1pripep_rename = SynGAP_workspace_Dir + '/renamed_id/' + str(sp1) + '.rename.primary.pep'
    sp2pripep_rename = SynGAP_workspace_Dir + '/renamed_id/' + str(sp2) + '.rename.primary.pep'
    os.system('python -m jcvi.formats.gff bed --type=' + str(annoType1) +
              ' --key=' +str(annoKey1) + ' --parent_key=' + str(annoparentKey1) +
              ' --primary_only ' + sp1gff_rename + ' -o ' + sp1bed_rename)
    os.system('python -m jcvi.formats.gff bed --type=' + str(annoType2) +
              ' --key=' + str(annoKey2) + ' --parent_key=' + str(annoparentKey2) +
              ' --primary_only ' + sp2gff_rename + ' -o ' + sp2bed_rename)
    os.system('cut -f4 ' + sp1bed_rename + ' > ' + sp1priid_rename)
    os.system('cut -f4 ' + sp2bed_rename + ' > ' + sp2priid_rename)
    os.system('gffread ' + sp1gff_rename + ' -g ' + sp1fa + ' -x ' + sp1cds_rename)
    os.system('seqkit grep -f ' + sp1priid_rename + ' ' + sp1cds_rename + ' > ' + sp1pricds_rename)
    os.system('gffread ' + sp1gff_rename + ' -g ' + sp1fa + ' -y ' + sp1pep_rename + ' -S')
    os.system('seqkit grep -f ' + sp1priid_rename + ' ' + sp1pep_rename + ' > ' + sp1pripep_rename)
    os.system('gffread ' + sp2gff_rename + ' -g ' + sp2fa + ' -x ' + sp2cds_rename)
    os.system('seqkit grep -f ' + sp2priid_rename + ' ' + sp2cds_rename + ' > ' + sp2pricds_rename)
    os.system('gffread ' + sp2gff_rename + ' -g ' + sp2fa + ' -y ' + sp2pep_rename + ' -S')
    os.system('seqkit grep -f ' + sp2priid_rename + ' ' + sp2pep_rename + ' > ' + sp2pripep_rename)
    ## perform MCScan using renamed_id files
    global mcscan_type
    if str(datatype) == 'nucl':
        mcscan_type = 'mcscan_cds'
        os.mkdir(SynGAP_workspace_Dir + '/' + mcscan_type)
        os.system('ln -s ' + SynGAP_workspace_Dir + '/renamed_id/' + str(sp1) + '.rename.primary.cds ' +
                      SynGAP_workspace_Dir + '/' + mcscan_type + '/' + str(sp1) + '.cds')
        os.system('ln -s ' + SynGAP_workspace_Dir + '/renamed_id/' + str(sp1) + '.rename.bed ' +
                      SynGAP_workspace_Dir + '/' + mcscan_type + '/' + str(sp1) + '.bed')
        os.system('ln -s ' + SynGAP_workspace_Dir + '/renamed_id/' + str(sp2) + '.rename.primary.cds ' +
                      SynGAP_workspace_Dir + '/' + mcscan_type + '/' + str(sp2) + '.cds')
        os.system('ln -s ' + SynGAP_workspace_Dir + '/renamed_id/' + str(sp2) + '.rename.bed ' +
                      SynGAP_workspace_Dir + '/' + mcscan_type + '/' + str(sp2) + '.bed')
    elif str(datatype) == 'prot':
        mcscan_type = 'mcscan_pep'
        os.mkdir(SynGAP_workspace_Dir + '/' + mcscan_type)
        os.system('ln -s ' + SynGAP_workspace_Dir + '/renamed_id/' + str(sp1) + '.rename.primary.pep ' +
                      SynGAP_workspace_Dir + '/' + mcscan_type + '/' + str(sp1) + '.pep')
        os.system('ln -s ' + SynGAP_workspace_Dir + '/renamed_id/' + str(sp1) + '.rename.bed ' +
                      SynGAP_workspace_Dir + '/' + mcscan_type + '/' + str(sp1) + '.bed')
        os.system('ln -s ' + SynGAP_workspace_Dir + '/renamed_id/' + str(sp2) + '.rename.primary.pep ' +
                      SynGAP_workspace_Dir + '/' + mcscan_type + '/' + str(sp2) + '.pep')
        os.system('ln -s ' + SynGAP_workspace_Dir + '/renamed_id/' + str(sp2) + '.rename.bed ' +
                      SynGAP_workspace_Dir + '/' + mcscan_type + '/' + str(sp2) + '.bed')
    os.chdir(SynGAP_workspace_Dir + '/' + mcscan_type)
    os.system('python -m jcvi.compara.catalog ortholog ' + str(sp1) + ' ' + str(sp2) +
                     ' --dbtype=' + str(datatype) + ' --cscore='  + str(cscore) + ' --cpus=' +
                     str(threads) + ' --no_strip_names --notex')
    print('\n[\033[0;36mINFO\033[0m]  MCScan for \033[0;35m' + sp1 + '\033[0m '
          'and \033[0;35m' + sp2 + '\033[0m Done!\n')

    # Find gaps from the .anchors file
    os.chdir(workDir)
    print('[\033[0;36mINFO\033[0m] Finding gaps from the .anchors file, please wait ...')
    anchors = SynGAP_workspace_Dir + '/' + str(sp1) + '.' + str(sp2) + '.anchors'
    anchors_gap = SynGAP_workspace_Dir + '/' + str(sp1) + '.' + str(sp2) + '.anchors.gap'
    os.system('ln -s ' + SynGAP_workspace_Dir + '/' + mcscan_type + '/' + str(
        sp1) + '.' + str(sp2) + '.anchors ' + anchors)
    print('[\033[0;36mINFO\033[0m] Finding ...')
    gapfind.gap_find(sp1bed_rename, sp2bed_rename, anchors, anchors_gap)
    print('[\033[0;36mINFO\033[0m] Running Done!\n')

    # Name the renamed geneids back to origin geneids
    print('[\033[0;36mINFO\033[0m] Naming the renamed geneids back to original geneids in '
          '`\033[0;35m.anchors\033[0m` and `\033[0;35m.anchors.gap\033[0m`, please wait ...')
    anchors_originalid = SynGAP_workspace_Dir + '/' + str(sp1) + '.' + str(
        sp2) + '.originalid.anchors'
    anchors_gap_originalid = SynGAP_workspace_Dir + '/' + str(sp1) + '.' + str(
        sp2) + '.originalid.anchors.gap'
    nameback.nameback(anchors, sp1mRNAmap, sp2mRNAmap, anchors_originalid)
    print('[\033[0;36mINFO\033[0m] Running Done!')
    nameback.nameback(anchors_gap, sp1mRNAmap, sp2mRNAmap, anchors_gap_originalid)
    print('[\033[0;36mINFO\033[0m] Running Done!\n')
    Sgenepairs_originalid = SynGAP_workspace_Dir + '/' + str(sp1) + '.' + str(sp2) + '.originalid.anchors.genepairs'
    print('[\033[0;36mINFO\033[0m] Extracting syntenic gene pairs from MCScan results for '
          '\033[0;35m' + str(sp1) + '\033[0m and \033[0;35m' + str(sp2) + '\033[0m , please wait ...')
    getgenepairs.get_genepair(anchors_originalid, Sgenepairs_originalid)
    print('[\033[0;36mINFO\033[0m] Running Done!\n')
    os.system('ln -s ' + anchors_originalid + ' ' + SynGAP_results_Dir + '/' + str(sp1) + '.' + str(
        sp2) + '.anchors')
    os.system('ln -s ' + anchors_gap_originalid + ' ' + SynGAP_results_Dir + '/' + str(sp1) + '.' + str(
        sp2) + '.anchors.gap')
    os.system('ln -s ' + Sgenepairs_originalid + ' ' + SynGAP_results_Dir + '/' + str(sp1) + '.' + str(
        sp2) + '.anchors.genepairs')

    # Extract .bed information for the gaps from .anchors.gap
    print('[\033[0;36mINFO\033[0m] Extracting .bed information for the gaps, please wait ...')
    sp1gapbed = SynGAP_workspace_Dir + '/' + str(sp1) + '.gap.bed'
    sp2gapbed = SynGAP_workspace_Dir + '/' + str(sp2) + '.gap.bed'
    gapbed.gapbed(anchors_gap, sp1bed_rename, sp2bed_rename, sp1gapbed, sp2gapbed)
    print('[\033[0;36mINFO\033[0m] Running Done!\n')

    # Prepare squences data for gap annotation
    print('[\033[0;36mINFO\033[0m] Preparing sequences for gap annotation, please wait ...')
    os.mkdir(SynGAP_workspace_Dir + '/gapanno')
    os.chdir(SynGAP_workspace_Dir + '/gapanno')
    gapseq.gapseq(sp1fa, sp2fa, sp1pep_rename, sp2pep_rename, sp1gapbed, sp2gapbed)
    print('[\033[0;36mINFO\033[0m] Running Done!\n')

    # Perform gap annotation
    print('[\033[0;36mINFO\033[0m] Performing gap annotation, please wait ...')
    os.system('ln -s ' + str(script_dir) + '/bin/genBlast_v138_linux_x86_64/formatdb')
    os.system('ln -s ' + str(script_dir) + '/bin/genBlast_v138_linux_x86_64/blastall')
    os.system('ln -s ' + str(script_dir) + '/bin/genBlast_v138_linux_x86_64/alignscore.txt')
    gapanno_genblastg.gapanno(str(script_dir) + '/bin/genBlast_v138_linux_x86_64', str(evalue), str(rank),
                              str(coverage), str(threads))
    print('[\033[0;36mINFO\033[0m] Running Done!\n')

    # Split the modified annotation file
    print('[\033[0;36mINFO\033[0m] Spliting the raw modified annotations file, please wait ...')
    os.chdir(workDir)
    modgff = SynGAP_workspace_Dir + '/' + 'modified.raw.gff'
    sp1modgff = SynGAP_workspace_Dir + '/' + str(sp1) + '.modified.raw.gff'
    sp2modgff = SynGAP_workspace_Dir + '/' + str(sp2) + '.modified.raw.gff'
    splitmodgff_genblastg.split(str(sp1), str(sp2), sp1gapbed, sp2gapbed, modgff, sp1modgff, sp2modgff)
    print('[\033[0;36mINFO\033[0m] Running Done!\n')

    # Name the renamed geneids back to original geneids in modified annotation files
    print('[\033[0;36mINFO\033[0m] Naming the renamed geneids back to '
          'original geneids in modified annotation files, please wait ...')
    sp1modgff_nameback = SynGAP_workspace_Dir + '/' + str(sp1) + '.modified.gff'
    sp2modgff_nameback = SynGAP_workspace_Dir + '/' + str(sp2) + '.modified.gff'
    gffnameback.nameback_gff(sp1modgff, sp2mRNAmap, sp1modgff_nameback)
    print('[\033[0;36mINFO\033[0m] Running Done!')
    gffnameback.nameback_gff(sp2modgff, sp1mRNAmap, sp2modgff_nameback)
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
        'python -m jcvi.formats.gff bed --type=mRNA --key=ID --parent_key=Parent ' + sp1modgff_nameback + ' -o ' + sp1modbed)
    os.system(
        'python -m jcvi.formats.gff bed --type=mRNA --key=ID --parent_key=Parent ' + sp2modgff_nameback + ' -o ' + sp2modbed)
    redundantfilter.redundantfilter(sp1modgff_nameback, sp1modbed, sp1modgff_disredundant)
    print('[\033[0;36mINFO\033[0m] Running Done!')
    redundantfilter.redundantfilter(sp2modgff_nameback, sp2modbed, sp2modgff_disredundant)
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
    os.system("grep mRNA " + sp1modgff_filtered + " | cut -f9 | cut -d';' -f2 | sed 's/Name=" + str(sp1) +
              "_//g' | sed 's/_[1-9][0-9]$//g' | sed 's/_[1-9]$//g' | sort | uniq > " + sp1modpep_originalid)
    os.system("grep mRNA " + sp2modgff_filtered + " | cut -f9 | cut -d';' -f2 | sed 's/Name=" + str(sp2) +
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
    Rlabel.Rlabel(str(sp1), sp1modgff_filtered, sp1R, sp1modgff_R)
    print('[\033[0;36mINFO\033[0m] Running Done!')
    print('[\033[0;36mINFO\033[0m] Calculating and labeling \033[0;33mR value\033[0m '
          'for modified annotations of \033[0;35m' + str(sp2) + '\033[0m, please wait ...')
    Rlabel.Rlabel(str(sp2), sp2modgff_filtered, sp2R, sp2modgff_R)
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
    print('[\033[0;36mINFO\033[0m] Performing gene function prediction for original genome annotations, please wait ...')
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
    sp1modgff_miss_annotated_I = SynGAP_workspace_Dir + '/cluster/' + str(sp1) + '.modified.filtered.R.miss_annotated.I.gff'
    sp1modgff_mis_annotated_I = SynGAP_workspace_Dir + '/cluster/' + str(sp1) + '.modified.filtered.R.mis_annotated.I.gff'
    sp1modgff_miss_annotated_III = SynGAP_workspace_Dir + '/cluster/' + str(sp1) + '.modified.filtered.R.miss_annotated.III.gff'
    sp1modgff_mis_annotated_III = SynGAP_workspace_Dir + '/cluster/' + str(sp1) + '.modified.filtered.R.mis_annotated.III.gff'
    sp1modgff_miss_annotated_IIa = SynGAP_workspace_Dir + '/cluster/' + str(sp1) + '.modified.filtered.R.miss_annotated.IIa.gff'
    sp1modgff_mis_annotated_IIa = SynGAP_workspace_Dir + '/cluster/' + str(sp1) + '.modified.filtered.R.mis_annotated.IIa.gff'
    sp1modgff_miss_annotated_IIb = SynGAP_workspace_Dir + '/cluster/' + str(sp1) + '.modified.filtered.R.miss_annotated.IIb.gff'
    sp1modgff_mis_annotated_IIb = SynGAP_workspace_Dir + '/cluster/' + str(sp1) + '.modified.filtered.R.mis_annotated.IIb.gff'
    sp2modgff_miss_annotated_I = SynGAP_workspace_Dir + '/cluster/' + str(sp2) + '.modified.filtered.R.miss_annotated.I.gff'
    sp2modgff_mis_annotated_I = SynGAP_workspace_Dir + '/cluster/' + str(sp2) + '.modified.filtered.R.mis_annotated.I.gff'
    sp2modgff_miss_annotated_III = SynGAP_workspace_Dir + '/cluster/' + str(sp2) + '.modified.filtered.R.miss_annotated.III.gff'
    sp2modgff_mis_annotated_III = SynGAP_workspace_Dir + '/cluster/' + str(sp2) + '.modified.filtered.R.mis_annotated.III.gff'
    sp2modgff_miss_annotated_IIa = SynGAP_workspace_Dir + '/cluster/' + str(sp2) + '.modified.filtered.R.miss_annotated.IIa.gff'
    sp2modgff_mis_annotated_IIa = SynGAP_workspace_Dir + '/cluster/' + str(sp2) + '.modified.filtered.R.mis_annotated.IIa.gff'
    sp2modgff_miss_annotated_IIb = SynGAP_workspace_Dir + '/cluster/' + str(sp2) + '.modified.filtered.R.miss_annotated.IIb.gff'
    sp2modgff_mis_annotated_IIb = SynGAP_workspace_Dir + '/cluster/' + str(sp2) + '.modified.filtered.R.mis_annotated.IIb.gff'
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
    os.system('ln -s ' + sp1modgff_I + ' ' + sp1GAPcleangff)
    os.system('ln -s ' + sp2modgff_I + ' ' + sp2GAPcleangff)
    os.system('ln -s ' + sp1modgff_miss_annotated_I + ' ' + sp1GAPcleanmiss_annotatedgff)
    os.system('ln -s ' + sp1modgff_mis_annotated_I + ' ' + sp1GAPcleanmis_annotatedgff)
    os.system('ln -s ' + sp2modgff_miss_annotated_I + ' ' + sp2GAPcleanmiss_annotatedgff)
    os.system('ln -s ' + sp2modgff_mis_annotated_I + ' ' + sp2GAPcleanmis_annotatedgff)
    os.system('cat ' + sp1gff + ' ' + sp1modgff_I + ' > ' + sp1GAPgff)
    os.system('cat ' + sp2gff + ' ' + sp2modgff_I + ' > ' + sp2GAPgff)
    print('[\033[0;36mINFO\033[0m] Running Done!\n')

    # Perform MCScan for species1 and species2 after SynGAP
    print('[\033[0;36mINFO\033[0m] Performing MCScan for \033[0;35m' + str(sp1) + '\033[0m '
          'and \033[0;35m' + str(sp2) + '\033[0m after SynGAP, please wait ...')
    print('[\033[0;36mINFO\033[0m] Preparing files ...')
    mcscanafterSynGAP_Dir = SynGAP_workspace_Dir + '/mcscanafterSynGAP'
    os.mkdir(mcscanafterSynGAP_Dir)
    sp1GAPbed = mcscanafterSynGAP_Dir + '/' + str(sp1) + '.SynGAP.bed'
    sp2GAPbed = mcscanafterSynGAP_Dir + '/' + str(sp2) + '.SynGAP.bed'
    sp1GAPpriid = mcscanafterSynGAP_Dir + '/' + str(sp1) + '.SynGAP.primary.id'
    sp2GAPpriid = mcscanafterSynGAP_Dir + '/' + str(sp2) + '.SynGAP.primary.id'
    sp1GAPcds = mcscanafterSynGAP_Dir + '/' + str(sp1) + '.SynGAP.cds'
    sp2GAPcds = mcscanafterSynGAP_Dir + '/' + str(sp2) + '.SynGAP.cds'
    sp1GAPpep = mcscanafterSynGAP_Dir + '/' + str(sp1) + '.SynGAP.pep'
    sp2GAPpep = mcscanafterSynGAP_Dir + '/' + str(sp2) + '.SynGAP.pep'
    sp1GAPpricds = mcscanafterSynGAP_Dir + '/' + str(sp1) + '.SynGAP.primary.cds'
    sp2GAPpricds = mcscanafterSynGAP_Dir + '/' + str(sp2) + '.SynGAP.primary.cds'
    sp1GAPpripep = mcscanafterSynGAP_Dir + '/' + str(sp1) + '.SynGAP.primary.pep'
    sp2GAPpripep = mcscanafterSynGAP_Dir + '/' + str(sp2) + '.SynGAP.primary.pep'
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
    if str(datatype) == 'nucl':
        mcscan_type = 'mcscan_cds'
        os.mkdir(mcscanafterSynGAP_Dir + '/' + mcscan_type)
        os.system(
            'ln -s ' + sp1GAPpricds + ' ' + mcscanafterSynGAP_Dir + '/' + mcscan_type + '/' + str(sp1) + '.cds')
        os.system(
            'ln -s ' + sp1GAPbed + ' ' + mcscanafterSynGAP_Dir + '/' + mcscan_type + '/' + str(sp1) + '.bed')
        os.system(
            'ln -s ' + sp2GAPpricds + ' ' + mcscanafterSynGAP_Dir + '/' + mcscan_type + '/' + str(sp2) + '.cds')
        os.system(
            'ln -s ' + sp2GAPbed + ' ' + mcscanafterSynGAP_Dir + '/' + mcscan_type + '/' + str(sp2) + '.bed')
    elif str(datatype) == 'prot':
        mcscan_type = 'mcscan_pep'
        os.mkdir(mcscanafterSynGAP_Dir + '/' + mcscan_type)
        os.system(
            'ln -s ' + sp1GAPpripep + ' ' + mcscanafterSynGAP_Dir + '/' + mcscan_type + '/' + str(sp1) + '.pep')
        os.system(
            'ln -s ' + sp1GAPbed + ' ' + mcscanafterSynGAP_Dir + '/' + mcscan_type + '/' + str(sp1) + '.bed')
        os.system(
            'ln -s ' + sp2GAPpripep + ' ' + mcscanafterSynGAP_Dir + '/' + mcscan_type + '/' + str(sp2) + '.pep')
        os.system(
            'ln -s ' + sp2GAPbed + ' ' + mcscanafterSynGAP_Dir + '/' + mcscan_type + '/' + str(sp2) + '.bed')
    os.chdir(mcscanafterSynGAP_Dir + '/' + mcscan_type)
    os.system('python -m jcvi.compara.catalog ortholog ' + str(sp1) + ' ' + str(sp2) +
                     ' --dbtype=' + str(datatype) + ' --cscore='  + str(cscore) + ' --cpus=' +
                     str(threads) + ' --no_strip_names --notex')
    print('\n[\033[0;36mINFO\033[0m] MCScan for \033[0;35m' + str(sp1) + '\033[0m '
          'and \033[0;35m' + str(sp2) + '\033[0m after SynGAP Done!\n')

    # Extract syntenic gene pairs from MCScan results for species1 and species2 after SynGAP
    anchors_GAP = mcscanafterSynGAP_Dir + '/' + mcscan_type + '/' + str(
        sp1) + '.' + str(sp2) + '.anchors'
    Sgenepairs_GAP = mcscanafterSynGAP_Dir + '/' + mcscan_type + '/' + str(
        sp1) + '.' + str(sp2) + '.anchors.SynGAP.genepairs'
    print('[\033[0;36mINFO\033[0m] Extracting syntenic gene pairs from MCScan results for '
          '\033[0;35m' + str(sp1) + '\033[0m and \033[0;35m' + str(sp2) + '\033[0m '
          'after SynGAP, please wait ...')
    getgenepairs.get_genepair(anchors_GAP, Sgenepairs_GAP)
    os.chdir(SynGAP_results_Dir)
    os.system('ln -s ' + Sgenepairs_GAP)
    print('[\033[0;36mINFO\033[0m] Running Done!\n')

    # Finished message
    print('[\033[0;36mINFO\033[0m] SynGAP analysis for \033[0;35m' + str(sp1) + '\033[0m '
          'and \033[0;35m' + str(sp2) +'\033[0m Done!')
    print('[\033[0;36mINFO\033[0m] Please check the result files in `\033[0;35m' + SynGAP_results_Dir + '\033[0m`\n')
    os.chdir(workDir)
    return(SynGAP_results_Dir)


if __name__ == '__main__':
    sp1fa = sys.argv[1]
    sp1gff = sys.argv[2]
    sp2fa = sys.argv[3]
    sp2gff  = sys.argv[4]
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
    evalue = sys.argv[16]
    rank = sys.argv[17]
    coverage = sys.argv[18]
    SynGAP(sp1, sp2, annoType1, annoKey1, annoparentKey1, annoType2, annoKey2, annoparentKey2,
           datatype, cscore, threads, evalue, rank, coverage)