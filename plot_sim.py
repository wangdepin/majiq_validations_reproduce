import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
import argparse
import csv
from itertools import islice
from tools import operator, all_stats
import sys
import os

__author__    = "Jorge D. Vaquero"
__copyright__ = "Licensed under GNU GPLv3"


valtypes = ['TOTAL_CHG', 'TOTAL_NO_CHG', 'FDR', 'FNR', 'MCC']
lims = {'FDR': ([0, 0.012], [0.10, 0.35]),
        'FNR': ([0, 0.10], [0.3, 0.8])}
def parse_inputfile(filelst, tools):


    res = {xx:{} for xx in valtypes}
    res_std = {xx: {} for xx in valtypes}
    nfls = len(filelst)
    tlist = set()
    for fidx, fname in enumerate(filelst):
        fp = islice(open(fname, "r"), 1, None)
        for row in csv.DictReader(fp, delimiter='\t'):
            if row['Tools'] not in tools:
                continue
            for xx in valtypes:
                v, vstd = row[xx].split('(')
                vstd = float(vstd.split(')')[0])
                v = float(v)
                try:
                    res[xx][row['Tools']][fidx] = v
                    res_std[xx][row['Tools']][fidx] = vstd
                except KeyError:
                    res[xx][row['Tools']] = np.zeros(nfls)
                    res[xx][row['Tools']][fidx] = v
                    res_std[xx][row['Tools']] = np.zeros(nfls)
                    res_std[xx][row['Tools']][fidx] = vstd
                tlist.add(row['Tools'])

    # print ('KK', res_std)
    return res, res_std, list(tlist)

def main(ifile, outfile, tools):
    if tools is None:
        tools = all_stats

    res, res_std, tlist = parse_inputfile(ifile, tools)
    print(res)
    clist = {}
    for xxi in tlist:
        if xxi == 'whippet_complex':
            xx = 'whippet'
        else:
            xx = xxi;
        module_ = __import__('tools.' + xx.lower(), fromlist=xx.title())
        class_ = getattr(module_, xx.title())
        #operator[xx] = class_()
        color, lst = class_().get_color()
        clist[xxi] = color

    width = 1
    global_space = len(tlist) * width + 0.4
    xticks = np.arange(global_space/2, (global_space+1)*len(ifile)-(global_space/2), global_space)


    print(res)
    print(len(res.keys()))
    num_subplots = len(res.keys())+2
    fig = plt.figure(figsize=(10, 10))
    gs = gridspec.GridSpec(num_subplots, 1, height_ratios=[2, 2, 1, 1, 1, 1, 2])

    #fig, ax = plt.subplots(num_subplots, 1, sharex=True, figsize=(10, 10))
    ax_idx = 0
    ax = [None]*num_subplots
    ax[-1] = fig.add_subplot(gs[num_subplots-1])
    for x in range(num_subplots-1):
        ax[x] = fig.add_subplot(gs[x], sharex=ax[-1])
        plt.setp(ax[x].get_xticklabels(), visible=False)

    for valtype, values in res.items():
        idx = 0
        for tid, vals in values.items():
            initx = idx * width
            xval = np.arange(initx, initx+(global_space*len(vals)), global_space)[:len(vals)]

            ax[ax_idx].bar(xval, vals, color=clist[tid], label=tid, width=1)
            if valtype in lims:
                ax[ax_idx+1].bar(xval, vals, color=clist[tid], label=tid, width=1)
                ax[ax_idx].set_ylim(lims[valtype][1])
                ax[ax_idx+1].set_ylim(lims[valtype][0])
                ax[ax_idx].spines['bottom'].set_visible(False)
                ax[ax_idx+1].spines['top'].set_visible(False)
                ax[ax_idx].xaxis.tick_top()
                ax[ax_idx].tick_params(labeltop=False)  # don't put tick labels at the top
                ax[ax_idx+1].xaxis.tick_bottom()

                for i, v in enumerate(vals):
                    if v > lims[valtype][0][1]:
                        ax[ax_idx+1].plot(xval[i], lims[valtype][0][1]-0.005, 'ro', color=clist[tid])

#            ax[ax_idx].errorbar(range(len(ifile)), vals, res_std[valtype][tid],linestyle='None', fmt='-o', label=tid, color=color)
            idx += 1
        ax[ax_idx].set_ylabel(valtype)
        ax_idx += 2 if valtype in lims else 1

#    ax[-1].set_xticklabels(ifile)
    plt.xticks(xticks)
    #print ("KK", ifile)
#    ax[-1].set_xticklabels([os.path.basename(xx) for xx in ifile])
    plt.legend()
    plt.savefig(outfile)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('input_file', nargs='+')
    parser.add_argument('-t', '--tools', nargs='+', default=None)
    parser.add_argument('-o', '--output', action='store', dest='outfile', default='./plot.pdf',
                        help='Output file name')
    # parser.add_argument('--labels', help='The labels for the plot lines of the ratios.')
    args = parser.parse_args()
    print(' '.join(sys.argv))

    filelist = args.input_file
    main(filelist, tools=args.tools, outfile=args.outfile)




