#!/usr/bin/env python
# -*- coding:utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from kneed import KneeLocator


def EVIplot(args):
    df = pd.read_csv(args.EVI, sep="\t", header=0)
    df['Rank'] = range(1, len(df['EVI']) + 1)
    ybg = np.array(df['EVI'])
    xbg = np.array(df["Rank"])

    y_max = max(df['EVI'])
    num_ticks = 4
    step = int((y_max - 0) / num_ticks)
    ticks = np.arange(0, y_max + step, step)

    fig_width = int(args.figsize.split('x')[0])
    fig_length = int(args.figsize.split('x')[1])
    plt.figure(figsize=(fig_width, fig_length))
    ax = plt.gca()
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    plt.xticks(fontsize=int(args.fontsize))
    plt.yticks(ticks, fontsize=int(args.fontsize))
    plt.scatter(xbg, ybg, s=int(args.fontsize) * 1.4, c="LightGrey")
    plt.xlabel('Rank', fontsize=int(args.fontsize) * 1.4, labelpad=int(args.fontsize))
    plt.ylabel('EVI', fontsize=int(args.fontsize) * 1.4, labelpad=int(args.fontsize))
    f1 = np.polyfit(xbg, ybg, 10)
    p1 = np.poly1d(f1)
    x_vals = list(range(1, 15000, 10))
    y_vals = p1(x_vals)
    kneedle = KneeLocator(x_vals, y_vals, S=1, curve="convex", direction="decreasing", online=True)
    plt.hlines(round(kneedle.knee_y, 2), 1, len(df['EVI']) + 1,
               colors='r', linestyles='dashed', label='EVI = ' + str(round(kneedle.knee_y, 2)))
    plt.legend(loc=1, frameon=False, fontsize=int(args.fontsize) * 1.4)
    if args.highlightid is not None:
        highlight_id = pd.read_csv(args.highlightid, sep="\t", header=0)
        for i in highlight_id['GeneID1']:
            hl_point = df.loc[df['GeneID1'] == i]
            y = hl_point['EVI']
            x = hl_point['Rank']
            label = str(highlight_id.loc[highlight_id['GeneID1'] == i]['Label']).strip().split()[1]
            plt.scatter(x, y, s=int(args.fontsize) * 1.4, c=args.highlightcolor)
            if label != 'NaN':
                plt.text(x + int(args.fontsize) * 30, y, label, fontsize=int(args.fontsize) * 1.4)
    plt.savefig(args.outgraph)
    plt.close()
