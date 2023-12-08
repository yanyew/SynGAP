#!/usr/bin/env python
# -*- coding:utf-8 -*-
import math
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from kneed import KneeLocator
from scipy.stats import pearsonr


def filter_exp(genelist, expfile):
    df = pd.read_csv(expfile, sep="\t", header=0, index_col="GeneID")
    df_new = df.reindex(set(genelist))
    df_new.to_csv(expfile + ".pair.tmp", sep='\t', header=True)


def EVI(args):
    genedict1 = {}
    genedict2 = {}
    geneset1 = []
    geneset2 = []
    gene_type = {}
    gene_type1 = {}
    count = 0
    with open(args.genepair) as genefile:
        for line in genefile.readlines():
            line = line.strip().split('\t')
            geneset1.append(line[0])
            geneset2.append(line[1])
            genedict1[count] = line[0]
            genedict2[count] = line[1]
            gene_type[line[0]] = line[2]
            if len(line) == 4:
                gene_type1[line[0]] = line[3]
            count += 1

    filter_exp(geneset1, args.sp1exp)
    filter_exp(geneset2, args.sp2exp)
    Expdataset1 = pd.read_csv(args.sp1exp + ".pair.tmp", sep="\t", header=0, index_col="GeneID")
    Expdataset2 = pd.read_csv(args.sp2exp + ".pair.tmp", sep="\t", header=0, index_col="GeneID")

    Expavgset1 = Expdataset1.mean(axis=1)
    Expavgset2 = Expdataset2.mean(axis=1)
    w1, w2, w3 = (args.weight).split(':')

    global output_tf, output_nontf
    if gene_type1:
        output = open(args.genepair + ".EVI.raw.xls", 'w+')
        output.write("GeneID1\tGeneID2\tHigherAvgExpGeneID\tGeneType\tAvgExpMax\tFoldchange\tPearson\tEVI\n")
        output_se = open(args.genepair + ".SEG.EVI.xls", 'w+')
        output_se.write("GeneID1\tGeneID2\tHigherAvgExpGeneID\tGeneType\tAvgExpMax\tFoldchange\tPearson\tEVI\n")
        output_tf = open(args.genepair + ".TF.EVI.xls", 'w+')
        output_tf.write("GeneID1\tGeneID2\tHigherAvgExpGeneID\tGeneType\tAvgExpMax\tFoldchange\tPearson\tEVI\n")
        output_nontf = open(args.genepair + ".nonTF.EVI.xls", 'w+')
        output_nontf.write("GeneID1\tGeneID2\tHigherAvgExpGeneID\tGeneType\tAvgExpMax\tFoldchange\tPearson\tEVI\n")
        for key in genedict1:
            geneid1 = genedict1[key]
            geneid2 = genedict2[key]
            higherexpgeneid = ''
            specificexpgeneid = []
            Foldchange = float()
            Expavgmax = max(Expavgset1.loc[geneid1], Expavgset2.loc[geneid2])
            Expavgmin = min(Expavgset1.loc[geneid1], Expavgset2.loc[geneid2])
            pearson = pearsonr(Expdataset1.loc[geneid1], Expdataset2.loc[geneid2])[0]
            warnings.simplefilter('ignore')
            if Expavgmax == Expavgset1.loc[geneid1]:
                higherexpgeneid = geneid1
                if Expavgmin < 0.1 and Expavgmax >= 0.1:
                    Foldchange = Expavgmax / 0.1
                    specificexpgeneid.append(geneid1)
                elif Expavgmin < 0.1 and Expavgmax < 0.1:
                    Foldchange = 1
                    pearson = 1
                else:
                    Foldchange = Expavgmax / Expavgmin
            elif Expavgmax == Expavgset2.loc[geneid2]:
                higherexpgeneid = geneid2
                if Expavgmin < 0.1 and Expavgmax >= 0.1:
                    Foldchange = Expavgmax / 0.1
                    specificexpgeneid.append(geneid1)
                elif Expavgmin < 0.1 and Expavgmax < 0.1:
                    Foldchange = 1
                    pearson = 1
                else:
                    Foldchange = Expavgmax / Expavgmin
            if pearson >= -1 and pearson <= 1:
                EVI = math.log2((Expavgmax + 1) ** int(w1) * Foldchange ** int(w2) * (2 - pearson) ** int(w3))
                output.write(
                    geneid1 + "\t" + geneid2 + "\t" + higherexpgeneid + "\t" + gene_type1[geneid1] + "\t" +
                    str(round(Expavgmax, 2)) + "\t" + str(round(Foldchange, 2)) + "\t" + \
                    str(round(pearson, 2)) + "\t" + str(round(EVI, 2)) + "\n")
                if geneid1 in specificexpgeneid:
                    output_se.write(
                        geneid1 + "\t" + geneid2 + "\t" + higherexpgeneid + "\t" + gene_type1[geneid1] + "\t" +
                        str(round(Expavgmax, 2)) + "\t" + str(round(Foldchange, 2)) + "\t" + \
                        str(round(pearson, 2)) + "\t" + str(round(EVI, 2)) + "\n")
                if gene_type1[geneid1] == 'nonTF/TR':
                    output_nontf.write(
                        geneid1 + "\t" + geneid2 + "\t" + higherexpgeneid + "\t" + gene_type1[geneid1] + "\t" +
                        str(round(Expavgmax, 2)) + "\t" + str(round(Foldchange, 2)) + "\t" + \
                        str(round(pearson, 2)) + "\t" + str(round(EVI, 2)) + "\n")
                elif gene_type1[geneid1] != 'nonTF/TR':
                    output_tf.write(
                        geneid1 + "\t" + geneid2 + "\t" + higherexpgeneid + "\t" + gene_type1[geneid1] + "\t" +
                        str(round(Expavgmax, 2)) + "\t" + str(round(Foldchange, 2)) + "\t" + \
                        str(round(pearson, 2)) + "\t" + str(round(EVI, 2)) + "\n")
            else:
                pearson1 = 0
                EVI = math.log2((Expavgmax + 1) ** int(w1) * Foldchange ** int(w2) * (2 - pearson1) ** int(w3))
                output.write(
                    geneid1 + "\t" + geneid2 + "\t" + higherexpgeneid + "\t" + gene_type1[geneid1] + "\t" +
                    str(round(Expavgmax, 2)) + "\t" + str(round(Foldchange, 2)) + "\t" + \
                    str(round(pearson1, 2)) + "\t" + str(round(EVI, 2)) + "\n")
                if geneid1 in specificexpgeneid:
                    output_se.write(
                        geneid1 + "\t" + geneid2 + "\t" + higherexpgeneid + "\t" + gene_type1[geneid1] + "\t" +
                        str(round(Expavgmax, 2)) + "\t" + str(round(Foldchange, 2)) + "\t" + \
                        str(round(pearson1, 2)) + "\t" + str(round(EVI, 2)) + "\n")
                if gene_type1[geneid1] == 'nonTF/TR':
                    output_nontf.write(
                        geneid1 + "\t" + geneid2 + "\t" + higherexpgeneid + "\t" + gene_type1[geneid1] + "\t" +
                        str(round(Expavgmax, 2)) + "\t" + str(round(Foldchange, 2)) + "\t" + \
                        str(round(pearson, 2)) + "\t" + str(round(EVI, 2)) + "\n")
                elif gene_type1[geneid1] != 'nonTF/TR':
                    output_tf.write(
                        geneid1 + "\t" + geneid2 + "\t" + higherexpgeneid + "\t" + gene_type1[geneid1] + "\t" +
                        str(round(Expavgmax, 2)) + "\t" + str(round(Foldchange, 2)) + "\t" + \
                        str(round(pearson, 2)) + "\t" + str(round(EVI, 2)) + "\n")
        output.close()
        output_tf.close()
        output_nontf.close()
    else:
        output = open(args.genepair + ".EVI.raw.xls", 'w+')
        output.write("GeneID1\tGeneID2\tHigherAvgExpGeneID\tGeneType\tAvgExpMax\tFoldchange\tPearson\tEVI\n")
        output_se = open(args.genepair + ".SEG.EVI.xls", 'w+')
        output_se.write("GeneID1\tGeneID2\tHigherAvgExpGeneID\tGeneType\tAvgExpMax\tFoldchange\tPearson\tEVI\n")
        for key in genedict1:
            geneid1 = genedict1[key]
            geneid2 = genedict2[key]
            higherexpgeneid = ''
            specificexpgeneid = []
            Foldchange = float()
            Expavgmax = max(Expavgset1.loc[geneid1], Expavgset2.loc[geneid2])
            Expavgmin = min(Expavgset1.loc[geneid1], Expavgset2.loc[geneid2])
            pearson = pearsonr(Expdataset1.loc[geneid1], Expdataset2.loc[geneid2])[0]
            warnings.simplefilter('ignore')
            if Expavgmax == Expavgset1.loc[geneid1]:
                higherexpgeneid = geneid1
                if Expavgmin < 0.1 and Expavgmax >= 0.1:
                    Foldchange = Expavgmax / 0.1
                    specificexpgeneid.append(geneid1)
                elif Expavgmin < 0.1 and Expavgmax < 0.1:
                    Foldchange = 1
                    pearson = 1
                else:
                    Foldchange = Expavgmax / Expavgmin
            elif Expavgmax == Expavgset2.loc[geneid2]:
                higherexpgeneid = geneid2
                if Expavgmin < 0.1 and Expavgmax >= 0.1:
                    Foldchange = Expavgmax / 0.1
                    specificexpgeneid.append(geneid1)
                elif Expavgmin < 0.1 and Expavgmax < 0.1:
                    Foldchange = 1
                    pearson = 1
                else:
                    Foldchange = Expavgmax / Expavgmin
            if pearson >= -1 and pearson <= 1:
                EVI = math.log2((Expavgmax + 1) ** int(w1) * Foldchange ** int(w2) * (2 - pearson) ** int(w3))
                output.write(
                    geneid1 + "\t" + geneid2 + "\t" + higherexpgeneid + "\t" + gene_type[geneid1] + "\t" +
                    str(round(Expavgmax, 2)) + "\t" + str(round(Foldchange, 2)) + "\t" + \
                    str(round(pearson, 2)) + "\t" + str(round(EVI, 2)) + "\n")
                if geneid1 in specificexpgeneid:
                    output_se.write(
                        geneid1 + "\t" + geneid2 + "\t" + higherexpgeneid + "\t" + gene_type[geneid1] + "\t" +
                        str(round(Expavgmax, 2)) + "\t" + str(round(Foldchange, 2)) + "\t" + \
                        str(round(pearson, 2)) + "\t" + str(round(EVI, 2)) + "\n")
            else:
                pearson1 = 0
                EVI = math.log2((Expavgmax + 1) ** int(w1) * Foldchange ** int(w2) * (2 - pearson1) ** int(w3))
                output.write(
                    geneid1 + "\t" + geneid2 + "\t" + higherexpgeneid + "\t" + gene_type[geneid1] + "\t" +
                    str(round(Expavgmax, 2)) + "\t" + str(round(Foldchange, 2)) + "\t" + \
                    str(round(pearson1, 2)) + "\t" + str(round(EVI, 2)) + "\n")
                if geneid1 in specificexpgeneid:
                    output_se.write(
                        geneid1 + "\t" + geneid2 + "\t" + higherexpgeneid + "\t" + gene_type[geneid1] + "\t" +
                        str(round(Expavgmax, 2)) + "\t" + str(round(Foldchange, 2)) + "\t" + \
                        str(round(pearson1, 2)) + "\t" + str(round(EVI, 2)) + "\n")
        output.close()

    dataset_all = pd.read_csv(args.genepair + ".EVI.raw.xls", sep="\t", header=0)
    df_all = dataset_all.sort_values(by='EVI', ascending=False)
    df1_all = df_all.reset_index(drop=True)
    df1_all.to_csv(args.genepair + ".EVI.raw.xls", sep='\t', index=False, header=True)
    df1_all = df1_all.drop_duplicates(subset=['GeneID1'], keep='last')
    df1_all = df1_all.reset_index(drop=True)
    df1_all.to_csv(args.genepair + ".EVI.xls", sep='\t', index=False, header=True)

    dataset_se = pd.read_csv(args.genepair + ".SEG.EVI.xls", sep="\t", header=0)
    df_se = dataset_se.sort_values(by='EVI', ascending=False)
    df1_se = df_se.reset_index(drop=True)
    df1_se = df1_se.drop_duplicates(subset=['GeneID1'], keep='last')
    df1_se = df1_se.reset_index(drop=True)
    df1_se.to_csv(args.genepair + ".SEG.EVI.xls", sep='\t', index=False, header=True)

    df1_all['Rank'] = range(1, len(df1_all['EVI']) + 1)
    x = np.array(df1_all['Rank'])
    y = np.array(df1_all['EVI'])
    plt.figure(figsize=(20, 10))
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.scatter(x, y, c="LightGrey")
    plt.xlabel('Rank', fontsize=20, labelpad=10)
    plt.ylabel('EVI', fontsize=20, labelpad=10)
    f1 = np.polyfit(x, y, 10)
    p1 = np.poly1d(f1)
    x_vals = list(range(1, 15000, 10))
    y_vals = p1(x_vals)
    kneedle = KneeLocator(x_vals, y_vals, S=1, curve="convex", direction="decreasing", online=True)
    output_threshold = open(args.genepair + ".EVI.threshold.txt", 'w+')
    output_threshold.write("Threshold for EVI: " + str(round(kneedle.knee_y, 2)))
    output_threshold.close()
    plt.hlines(round(kneedle.knee_y, 2), 1, len(df1_all['EVI']) + 1,
               colors='r', linestyles='dashed', label='EVI = ' + str(round(kneedle.knee_y, 2)))
    plt.legend(loc=1, frameon=False, fontsize=20)
    plt.savefig(args.genepair + ".EVI." + args.format)
    plt.close()

    df1_all['ML+1'] = df1_all['AvgExpMax'] + 1
    df1_all['(2-PCC)'] = 2 - df1_all['Pearson']
    df1_all['ML'] = np.log2(df1_all['ML+1']) * int(w1)
    df1_all['FC'] = np.log2(df1_all['Foldchange']) * int(w2)
    df1_all['PCC'] = np.log2(df1_all['(2-PCC)']) * int(w3)
    plt.figure(figsize=(20, 10))
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.bar(df1_all['Rank'], df1_all['ML'], label='log2((ML+1)^' + w1 + ')', color='#69CCE6')
    plt.bar(df1_all['Rank'], df1_all['FC'], bottom=df1_all['ML'], label='log2(FC^' + w2 + ')',
            color='#FF6666')
    plt.bar(df1_all['Rank'], df1_all['PCC'], bottom=df1_all['ML'] + df1_all['FC'],
            label='log2((2-PCC)^' + w3 + ')', color='#FFA870')
    plt.xlabel('Rank', fontsize=20, labelpad=10)
    plt.ylabel('EVI', fontsize=20, labelpad=10)
    plt.legend(loc=1, frameon=False, fontsize=20)
    plt.savefig(args.genepair + ".EVI.indexweight." + args.format)
    plt.close()

    df1_all['ML+1'] = df1_all['AvgExpMax'] + 1
    df1_all['(2-PCC)'] = 2 - df1_all['Pearson']
    df1_all['ML'] = np.log2(df1_all['ML+1']) * int(w1)
    df1_all['FC'] = np.log2(df1_all['Foldchange']) * int(w2)
    df1_all['PCC'] = np.log2(df1_all['(2-PCC)']) * int(w3)
    df1_all['EVIr'] = df1_all['ML'] + df1_all['FC'] + df1_all['PCC']
    df1_all['MLr'] = (df1_all['ML'] / df1_all['EVIr']) * 100
    df1_all['FCr'] = (df1_all['FC'] / df1_all['EVIr']) * 100
    df1_all['PCCr'] = (df1_all['PCC'] / df1_all['EVIr']) * 100
    plt.figure(figsize=(20, 10))
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.bar(df1_all['Rank'], df1_all['MLr'], label='log2((ML+1)^' + w1 + ')', color='#69CCE6')
    plt.bar(df1_all['Rank'], df1_all['FCr'], bottom=df1_all['MLr'], label='log2(FC^' + w2 + ')',
            color='#FF6666')
    plt.bar(df1_all['Rank'], df1_all['PCCr'], bottom=df1_all['MLr'] + df1_all['FCr'],
            label='log2((2-PCC)^' + w3 + ')', color='#FFA870')
    plt.xlabel('Rank', fontsize=20, labelpad=10)
    plt.ylabel('EVI', fontsize=20, labelpad=10)
    plt.legend(loc=1, frameon=True, fontsize=20)
    plt.savefig(args.genepair + ".EVI.indexweightratio." + args.format)
    plt.close()

    df1_se['Rank'] = range(1, len(df1_se['EVI']) + 1)
    x3 = np.array(df1_se['Rank'])
    y3 = np.array(df1_se['EVI'])
    plt.figure(figsize=(20, 10))
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.scatter(x3, y3, c="LightGrey")
    plt.xlabel('Rank', fontsize=20, labelpad=10)
    plt.ylabel('EVI', fontsize=20, labelpad=10)
    f4 = np.polyfit(x3, y3, 10)
    p4 = np.poly1d(f4)
    x_vals = list(range(1, 20000, 10))
    y_vals = p4(x_vals)
    kneedle = KneeLocator(x_vals, y_vals, S=1, curve="convex", direction="decreasing", online=True)
    output_threshold = open(args.genepair + ".SEG.EVI.threshold.txt", 'w+')
    output_threshold.write("Threshold for EVI: " + str(round(kneedle.knee_y, 3)))
    output_threshold.close()
    plt.hlines(round(kneedle.knee_y, 3), 1, len(df1_se['EVI']) + 1,
               colors='r', linestyles='dashed', label='EVI = ' + str(round(kneedle.knee_y, 3)))
    plt.legend(loc=1, frameon=False, fontsize=20)
    plt.savefig(args.genepair + ".SEG.EVI." + args.format)
    plt.close()

    if gene_type1:
        dataset_tf = pd.read_csv(args.genepair + ".TF.EVI.xls", sep="\t", header=0)
        df_tf = dataset_tf.sort_values(by='EVI', ascending=False)
        df1_tf = df_tf.reset_index(drop=True)
        df1_tf = df1_tf.drop_duplicates(subset=['GeneID1'], keep='last')
        df1_tf = df1_tf.reset_index(drop=True)
        df1_tf.to_csv(args.genepair + ".TF.EVI.xls", sep='\t', index=False, header=True)

        dataset_nontf = pd.read_csv(args.genepair + ".nonTF.EVI.xls", sep="\t", header=0)
        df_nontf = dataset_nontf.sort_values(by='EVI', ascending=False)
        df1_nontf = df_nontf.reset_index(drop=True)
        df1_nontf = df1_nontf.drop_duplicates(subset=['GeneID1'], keep='last')
        df1_nontf = df1_nontf.reset_index(drop=True)
        df1_nontf.to_csv(args.genepair + ".nonTF.EVI.xls", sep='\t', index=False, header=True)

        df1_tf['Rank'] = range(1, len(df1_tf['EVI']) + 1)
        x1 = np.array(df1_tf['Rank'])
        y1 = np.array(df1_tf['EVI'])
        plt.figure(figsize=(20, 10))
        plt.xticks(fontsize=16)
        plt.yticks(fontsize=16)
        plt.scatter(x1, y1, c="LightGrey")
        plt.xlabel('Rank', fontsize=20, labelpad=10)
        plt.ylabel('EVI', fontsize=20, labelpad=10)
        f2 = np.polyfit(x1, y1, 10)
        p2 = np.poly1d(f2)
        x_vals = list(range(1, 20000, 10))
        y_vals = p2(x_vals)
        kneedle = KneeLocator(x_vals, y_vals, S=1, curve="convex", direction="decreasing", online=True)
        output_threshold = open(args.genepair + ".TF.EVI.threshold.txt", 'w+')
        output_threshold.write("Threshold for EVI: " + str(round(kneedle.knee_y, 3)))
        output_threshold.close()
        plt.hlines(round(kneedle.knee_y, 3), 1, len(df1_tf['EVI']) + 1,
                   colors='r', linestyles='dashed', label='EVI = ' + str(round(kneedle.knee_y, 3)))
        plt.legend(loc=1, frameon=False, fontsize=20)
        plt.savefig(args.genepair + ".TF.EVI." + args.format)
        plt.close()

        df1_nontf['Rank'] = range(1, len(df1_nontf['EVI']) + 1)
        x2 = np.array(df1_nontf['Rank'])
        y2 = np.array(df1_nontf['EVI'])
        plt.figure(figsize=(20, 10))
        plt.xticks(fontsize=16)
        plt.yticks(fontsize=16)
        plt.scatter(x2, y2, c="LightGrey")
        plt.xlabel('Rank', fontsize=20, labelpad=10)
        plt.ylabel('EVI', fontsize=20, labelpad=10)
        f3 = np.polyfit(x2, y2, 10)
        p3 = np.poly1d(f3)
        x_vals = list(range(1, 20000, 10))
        y_vals = p3(x_vals)
        kneedle = KneeLocator(x_vals, y_vals, S=1, curve="convex", direction="decreasing", online=True)
        output_threshold = open(args.genepair + ".nonTF.EVI.threshold.txt", 'w+')
        output_threshold.write("Threshold for EVI: " + str(round(kneedle.knee_y, 3)))
        output_threshold.close()
        plt.hlines(round(kneedle.knee_y, 3), 1, len(df1_nontf['EVI']) + 1,
                   colors='r', linestyles='dashed', label='EVI = ' + str(round(kneedle.knee_y, 3)))
        plt.legend(loc=1, frameon=False, fontsize=20)
        plt.savefig(args.genepair + ".nonTF.EVI." + args.format)
        plt.close()
