#!/usr/bin/env python3
# -*- coding:utf-8 -*-
import os
import sys
import datetime


def genepair_in_anchors(infile, outfile):
    inf = open(infile, 'r')
    outf = open(outfile, 'w')
    while True:
        line = inf.readline()
        if not line: break
        if line[0] != '#':
            outf.write(line.split('\t')[0] + '\t' + line.split('\t')[1] + '\t' + 'Synteny' + '\n')
    inf.close()
    outf.close()


def genelist_bed(infile):
    inf_bed = open(infile, 'r')
    bed_genelist= []
    while True:
        line = inf_bed.readline()
        if not line: break
        gene_id = line.split('\t')[3]
        if gene_id not in bed_genelist:
            bed_genelist.append(gene_id)
    inf_bed.close()
    return bed_genelist


def gene_withoutSyn(infile_bed1, infile_bed2, genepair_in_anchors):
    genepair_anchors = open(genepair_in_anchors, 'r')
    genelist1 = genelist_bed(infile_bed1)
    genelist2 = genelist_bed(infile_bed2)
    genelist_anchors1 = []
    genelist_anchors2 = []
    while True:
        line = genepair_anchors.readline()
        if not line: break
        gene_id1 = line.split('\t')[0]
        gene_id2 = line.split('\t')[1]
        if gene_id1 not in genelist_anchors1:
            genelist_anchors1.append(gene_id1)
        if gene_id2 not in genelist_anchors2:
            genelist_anchors2.append(gene_id2)
    genelist_withoutSyn1 = list(set(genelist1).difference(set(genelist_anchors1)))
    genelist_withoutSyn2 = list(set(genelist2).difference(set(genelist_anchors2)))
    genepair_anchors.close()
    return genelist_withoutSyn1, genelist_withoutSyn2


def dict_fa(infile):
    inf = open(infile, 'r')
    seqs = {}
    seq_id = ''
    while True:
        line = inf.readline()
        if not line: break
        if line[0] == '>':
            seq_id = line[1:].strip()
            seqs[seq_id] = ''
        else:
            seqs[seq_id] += line.strip()
    inf.close()
    return seqs


def prepareseq(sp1, sp2, inseq_sp1, inseq_sp2, infile_bed1, infile_bed2, genepair_in_anchors):
    seqs_sp1 = dict_fa(inseq_sp1)
    seqs_sp2 = dict_fa(inseq_sp2)
    genelist_withoutSyn1, genelist_withoutSyn2 = gene_withoutSyn(infile_bed1, infile_bed2, genepair_in_anchors)
    query_sp1 = sp1 + '.query.fa'
    query_sp2 = sp2 + '.query.fa'
    outf_query_sp1 = open(query_sp1, 'a')
    outf_query_sp2 = open(query_sp2, 'a')
    for id in genelist_withoutSyn1:
        if id in seqs_sp1:
            outf_query_sp1.write('>' + id + '\n' + seqs_sp1[id] + '\n')
    for id in genelist_withoutSyn2:
        if id in seqs_sp2:
            outf_query_sp2.write('>' + id + '\n' + seqs_sp2[id] + '\n')
    outf_query_sp1.close()
    outf_query_sp2.close()


def twoway_blast(sp1, sp2, datatype, evalue, pnumber):
    sp12sp2_blast = sp1 + '2' + sp2 + '.blast'
    sp22sp1_blast = sp2 + '2' + sp1 + '.blast'
    query_sp1 = sp1 + '.query.fa'
    query_sp2 = sp2 + '.query.fa'
    if datatype == 'nucl':
        os.system('makeblastdb -in ' + query_sp1 + ' -dbtype nucl -out PDB_' + sp1 + ' > ' + sp1 + '.makeblastdb.log 2>&1')
        os.system('makeblastdb -in ' + query_sp2 + ' -dbtype nucl -out PDB_' + sp2 + ' > ' + sp2 + '.makeblastdb.log 2>&1')
        os.system('blastn -query ' + sp1 + '.query.fa -db PDB_' + sp2 + ' -out ' + sp12sp2_blast +
                  ' -outfmt 6 -evalue ' + str(evalue) + ' -max_target_seqs 1000 -num_threads ' + pnumber)
        os.system('blastn -query ' + sp2 + '.query.fa -db PDB_' + sp1 + ' -out ' + sp22sp1_blast +
                  ' -outfmt 6 -evalue ' + str(evalue) + ' -max_target_seqs 1000 -num_threads ' + pnumber)
    elif datatype == 'prot':
        os.system('makeblastdb -in ' + query_sp1 + ' -dbtype prot -out PDB_' + sp1 + ' > ' + sp1 + '.makeblastdb.log 2>&1')
        os.system('makeblastdb -in ' + query_sp2 + ' -dbtype prot -out PDB_' + sp2 + ' > ' + sp2 + '.makeblastdb.log 2>&1')
        os.system('blastp -query ' + sp1 + '.query.fa -db PDB_' + sp2 + ' -out ' + sp12sp2_blast +
                  ' -outfmt 6 -evalue ' + str(evalue) + ' -max_target_seqs 1000 -num_threads ' + pnumber)
        os.system('blastp -query ' + sp2 + '.query.fa -db PDB_' + sp1 + ' -out ' + sp22sp1_blast +
                  ' -outfmt 6 -evalue ' + str(evalue) + ' -max_target_seqs 1000 -num_threads ' + pnumber)


def blast_genepair(sp1, sp2):
    sp12sp2_blast = sp1 + '2' + sp2 + '.blast'
    sp22sp1_blast = sp2 + '2' + sp1 + '.blast'
    inf_sp12sp2_blast = open(sp12sp2_blast, 'r')
    inf_sp22sp1_blast = open(sp22sp1_blast, 'r')
    blast_dict_sp1 = {}
    blast_dict_sp2 = {}
    while True:
        line = inf_sp12sp2_blast.readline()
        if not line: break
        if line != '':
            query_id = line.split('\t')[0]
            subject_id = line.split('\t')[1]
            if query_id not in blast_dict_sp1.keys():
                blast_dict_sp1[query_id] = subject_id
    while True:
        line = inf_sp22sp1_blast.readline()
        if not line: break
        if line != '':
            query_id = line.split('\t')[0]
            subject_id = line.split('\t')[1]
            if query_id not in blast_dict_sp2.keys():
                blast_dict_sp2[query_id] = subject_id
    inf_sp12sp2_blast.close()
    inf_sp22sp1_blast.close()
    return blast_dict_sp1, blast_dict_sp2


def genepair(args):
    script_dir, filename = os.path.split(os.path.abspath(sys.argv[0]))
    global workDir
    workDir = os.getcwd()
    mission_time = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    genepair_Dir = workDir + '/genepair_' + mission_time
    genepair_results_Dir = workDir + '/genepair_' + mission_time + '/results'
    os.mkdir(genepair_Dir)
    os.mkdir(genepair_results_Dir)
    sp1fa = genepair_Dir + '/' + str(args.sp1) + '.fa'
    sp2fa = genepair_Dir + '/' + str(args.sp2) + '.fa'
    sp1gff = genepair_Dir + '/' + str(args.sp1) + '.gff3'
    sp2gff = genepair_Dir + '/' + str(args.sp2) + '.gff3'
    os.chdir(workDir)
    print('[\033[0;36mINFO\033[0m] Renaming the input files, please wait ...')
    print('[\033[0;36mINFO\033[0m] '
          'Rename `\033[0;35m' + str(args.sp1fa) + '\033[0m` into `\033[0;35m' + sp1fa + '\033[0m`')
    os.system('cp ' + str(args.sp1fa) + ' ' + sp1fa)
    print('[\033[0;36mINFO\033[0m] '
          'Rename `\033[0;35m' + str(args.sp2fa) + '\033[0m` into `\033[0;35m' + sp2fa + '\033[0m`')
    os.system('cp ' + str(args.sp2fa) + ' ' + sp2fa)
    print('[\033[0;36mINFO\033[0m] '
          'Rename `\033[0;35m' + str(args.sp1gff) + '\033[0m` into `\033[0;35m' + sp1gff + '\033[0m`')
    os.system('cp ' + str(args.sp1gff) + ' ' + sp1gff)
    print('[\033[0;36mINFO\033[0m] '
          'Rename `\033[0;35m' + str(args.sp2gff) + '\033[0m` into `\033[0;35m' + sp2gff + '\033[0m`')
    os.system('cp ' + str(args.sp2gff) + ' ' + sp2gff)
    print('[\033[0;36mINFO\033[0m] Running Done!\n')

    # Perform MCScan for species1 and species2
    print('[\033[0;36mINFO\033[0m] Performing MCScan for \033[0;35m' + args.sp1 + '\033[0m '
          'and \033[0;35m' + args.sp2 + '\033[0m, please wait ...')
    print('[\033[0;36mINFO\033[0m] Preparing files ...')
    ## prepare files for original_id
    sp1bed = genepair_Dir + '/' + str(args.sp1) + '.bed'
    sp2bed = genepair_Dir + '/' + str(args.sp2) + '.bed'
    sp1priid = genepair_Dir + '/' + str(args.sp1) + '.primary.id'
    sp2priid = genepair_Dir + '/' + str(args.sp2) + '.primary.id'
    sp1cds = genepair_Dir + '/' + str(args.sp1) + '.cds'
    sp2cds = genepair_Dir + '/' + str(args.sp2) + '.cds'
    sp1pep = genepair_Dir + '/' + str(args.sp1) + '.pep'
    sp2pep = genepair_Dir + '/' + str(args.sp2) + '.pep'
    sp1pricds = genepair_Dir + '/' + str(args.sp1) + '.primary.cds'
    sp2pricds = genepair_Dir + '/' + str(args.sp2) + '.primary.cds'
    sp1pripep = genepair_Dir + '/' + str(args.sp1) + '.primary.pep'
    sp2pripep = genepair_Dir + '/' + str(args.sp2) + '.primary.pep'
    os.system('python -m jcvi.formats.gff bed --type=' + str(args.annoType1) +
              ' --key=' + str(args.annoKey1) + ' --parent_key=' + str(args.annoparentKey1) +
              ' --primary_only ' + sp1gff + ' -o ' + sp1bed)
    os.system('python -m jcvi.formats.gff bed --type=' + str(args.annoType2) +
              ' --key=' + str(args.annoKey2) + ' --parent_key=' + str(args.annoparentKey2) +
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

    ## perform MCScan using renamed_id files
    global mcscan_type
    if str(args.datatype) == 'nucl':
        mcscan_type = 'mcscan_cds'
        os.mkdir(genepair_Dir + '/' + mcscan_type)
        os.system('ln -s ' + genepair_Dir + '/' + str(args.sp1) + '.primary.cds ' +
                  genepair_Dir + '/' + mcscan_type + '/' + str(args.sp1) + '.cds')
        os.system('ln -s ' + genepair_Dir + '/' + str(args.sp1) + '.bed ' +
                  genepair_Dir + '/' + mcscan_type + '/' + str(args.sp1) + '.bed')
        os.system('ln -s ' + genepair_Dir + '/' + str(args.sp2) + '.primary.cds ' +
                  genepair_Dir + '/' + mcscan_type + '/' + str(args.sp2) + '.cds')
        os.system('ln -s ' + genepair_Dir + '/' + str(args.sp2) + '.bed ' +
                  genepair_Dir + '/' + mcscan_type + '/' + str(args.sp2) + '.bed')
    elif str(args.datatype) == 'prot':
        mcscan_type = 'mcscan_pep'
        os.mkdir(genepair_Dir + '/' + mcscan_type)
        os.system('ln -s ' + genepair_Dir + '/' + str(args.sp1) + '.primary.pep ' +
                  genepair_Dir + '/' + mcscan_type + '/' + str(args.sp1) + '.pep')
        os.system('ln -s ' + genepair_Dir + '/' + str(args.sp1) + '.bed ' +
                  genepair_Dir + '/' + mcscan_type + '/' + str(args.sp1) + '.bed')
        os.system('ln -s ' + genepair_Dir + '/' + str(args.sp2) + '.primary.pep ' +
                  genepair_Dir + '/' + mcscan_type + '/' + str(args.sp2) + '.pep')
        os.system('ln -s ' + genepair_Dir + '/' + str(args.sp2) + '.bed ' +
                  genepair_Dir + '/' + mcscan_type + '/' + str(args.sp2) + '.bed')
    os.chdir(genepair_Dir + '/' + mcscan_type)
    os.system('python -m jcvi.compara.catalog ortholog ' + str(args.sp1) + ' ' + str(args.sp2) +
              ' --dbtype=' + str(args.datatype) + ' --cscore=' + str(args.cscore) + ' --cpus=' +
              str(args.threads) + ' --no_strip_names --notex')
    print('\n[\033[0;36mINFO\033[0m]  MCScan for \033[0;35m' + args.sp1 + '\033[0m '
          'and \033[0;35m' + args.sp2 + '\033[0m Done!\n')

    os.chdir(genepair_Dir)
    print('[\033[0;36mINFO\033[0m] Extracting syntenic gene pairs for \033[0;35m' + args.sp1 + '\033[0m '
          'and \033[0;35m' + args.sp2 + '\033[0m, please wait ...')
    anchors = genepair_Dir + '/' + mcscan_type + '/' + str(args.sp1) + '.' + str(args.sp2) + '.anchors'
    os.system('ln -s ' + anchors)
    pairwithSyn = args.sp1 + '.' + args.sp2 + '.Synteny.genepair'
    genepair_in_anchors(anchors, pairwithSyn)
    print('[\033[0;36mINFO\033[0m] Extracting Done!\n')

    print('[\033[0;36mINFO\033[0m] Extracting two-way best-blast gene pairs for \033[0;35m' + args.sp1 + '\033[0m '
          'and \033[0;35m' + args.sp2 + '\033[0m, please wait ...')
    if str(args.datatype) == 'nucl':
        prepareseq(str(args.sp1), str(args.sp2), sp1pricds, sp2pricds, sp1bed, sp2bed, pairwithSyn)
        twoway_blast(str(args.sp1), str(args.sp2), str(args.datatype), str(args.evalue), str(args.threads))
    elif str(args.datatype) == 'prot':
        prepareseq(str(args.sp1), str(args.sp2), sp1pripep, sp2pripep, sp1bed, sp2bed, pairwithSyn)
        twoway_blast(str(args.sp1), str(args.sp2), str(args.datatype), str(args.evalue), str(args.threads))
    blast_dict1, blast_dict2 = blast_genepair(str(args.sp1), str(args.sp2))
    pairwith2wayblast = args.sp1 + '.' + args.sp2 + '.2wayblast.genepair'
    out_pairwith2wayblast = open(pairwith2wayblast, 'w')
    for id1 in blast_dict1.keys():
        id2 = blast_dict1[id1]
        if id2 in blast_dict2 and id1 == blast_dict2[id2]:
            out_pairwith2wayblast.write(id1 + '\t' + id2 + '\t' + 'best two-way blast' + '\n')
    out_pairwith2wayblast.close()
    print('[\033[0;36mINFO\033[0m] Extracting Done!\n')

    print('[\033[0;36mINFO\033[0m] Merging syntenic and two-way best-blast gene pairs for \033[0;35m' + args.sp1 + '\033[0m '
          'and \033[0;35m' + args.sp2 + '\033[0m, please wait ...')
    genepair = args.sp1 + '.' + args.sp2 + '.final.genepair'
    if args.iTAK == 'no':
        os.system('cat ' + pairwithSyn + ' ' + pairwith2wayblast + ' > ' + genepair)
        print('[\033[0;36mINFO\033[0m] Merging Done!\n')
    elif args.iTAK == 'yes':
        infa = open(sp1pripep, 'r')
        os.system('perl ' + str(script_dir) + '/bin/iTAK/iTAK.pl -p ' + args.threads + ' -o ' + args.sp1 + '_iTAK ' + sp1pricds)
        inf_tf_c = open(args.sp1 + '_iTAK/tf_classification.txt', 'r')
        gene_type = {}
        while True:
            line = inf_tf_c.readline()
            if not line: break
            if line != '':
                lines = line.rstrip('\n').split('\t')
                gene_type[lines[0]] = lines[2] + '->' + lines[3]
        gene_list = []
        while True:
            line = infa.readline()
            if not line: break
            if line[0] == '>':
                geneid = line.rstrip('\n')[1:]
                gene_list.append(geneid)
        for geneid in gene_list:
            if geneid not in gene_type:
                gene_type[geneid] = 'nonTF/TR'
        infa.close()
        inf_tf_c.close()
        inf_pairwithSyn = open(pairwithSyn, 'r')
        inf_pairwith2wayblast = open(pairwith2wayblast, 'r')
        outf_genepair = open(genepair, 'w')
        while True:
            line = inf_pairwithSyn.readline()
            if not line: break
            lines = line.rstrip('\n').split('\t')
            outf_genepair.write(lines[0] + '\t' + lines[1] + '\t' + lines[2] + '\t' + gene_type[lines[0]] + '\n')
        while True:
            line = inf_pairwith2wayblast.readline()
            if not line: break
            lines = line.rstrip('\n').split('\t')
            outf_genepair.write(lines[0] + '\t' + lines[1] + '\t' + lines[2] + '\t' + gene_type[lines[0]] + '\n')
        outf_genepair.close()
        print('[\033[0;36mINFO\033[0m] Merging Done!\n')

    result_pairwithSyn = genepair_Dir + '/' + pairwithSyn
    result_pairwith2wayblast = genepair_Dir + '/' + pairwith2wayblast
    result_genepair = genepair_Dir + '/' + genepair
    os.chdir(genepair_results_Dir)
    os.system('ln -s ' + result_pairwithSyn)
    os.system('ln -s ' + result_pairwith2wayblast)
    os.system('ln -s ' + result_genepair)
    print('[\033[0;36mINFO\033[0m] Please check the result files in `\033[0;35m' + genepair_results_Dir + '\033[0m`\n')