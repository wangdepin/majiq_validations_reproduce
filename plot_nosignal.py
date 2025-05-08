import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import argparse
from tools import operator, all_stats
import os
from statistics import median, pstdev

__author__    = "Jorge D. Vaquero"
__copyright__ = "Licensed under GNU GPLv3"

plt.rc('font', size=18)


def main(tools_id, tools, label, runs_list, outfile):


    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(20, 10))
    nmethods = len(tools_id)

    signal   = np.zeros(shape=nmethods)
    nosignal = np.zeros(shape=nmethods)

    width = 1
    global_space = len(tools_id) * width + 0.4
    print (runs_list)
    xticks = np.arange(global_space/2, (global_space)*len(runs_list)-(global_space/2), global_space)
    print(xticks)
    for idx, (tool, runid) in enumerate(tools_id):
        val1, val2 = tools[(tool, runid)]
        color, lst = operator[tool].get_color()
        vals = []
        std_vals = []
        for xx in runs_list :
            try:
                v = float(median(val1[xx]))
                st = pstdev(val1[xx])
            except KeyError:
                v = float(-1)
                st = float(0)

            vals.append(v)
            std_vals.append(st)

        initx = idx*width
        xval = np.arange(initx, initx+(global_space*len(vals)), global_space)[:len(vals)]
        print (xval, vals)
        ax1.bar(xval, vals, yerr=std_vals, label=runid, color=color, edgecolor = "none", width=width)
        for ii in range(len(vals)):
#            print (vals[ii])
            if vals[ii] >= 0 :
                ax1.text(xval[ii] - (width/2), min(450, vals[ii]+20), "%d" % vals[ii], fontsize=10, rotation=45)

        vals_s = []
        std_vals_s = []
        for xi, xx in enumerate(runs_list) :
            try:
                vlist = [float(i1) / i2 for i1, i2 in zip(val1[xx], val2[xx])]
#                vlist = val1[xx] / val2[xx]
                v = median(vlist)
                st = pstdev(vlist)
                # v = vals[xi]/float(median(val2[xx]))
                #if (vals[xi] == -1 and float(val2[xx]) == -1):
                #    v = -1
            except KeyError:
                v = float(-0.1)
                st = float(0)
            vals_s.append(v)
            std_vals_s.append(st)

#        vals_s = [vals[xi]/float(xx) for xi, xx in enumerate(val2)]
        ax2.bar(xval, vals_s, yerr=std_vals_s, label=runid, color=color, edgecolor = "none", width=width)
        for ii in range(len(vals_s)):
#            print vals_s[ii]
            v = vals_s[ii] * 100
            if v >= 0 :
                ax2.text(xval[ii] - (width/2), min(vals_s[ii] + 0.01, 0.25), "%.1f%%" % v, fontsize=10, rotation=45)

    ax1.set_ylim([0, 500])
    ax1.set_ylabel('N Nosignal')
    ax2.set_ylabel('Estimated IIR')
    ax2.set_ylim([0, 0.3])
    ax2.set_xticklabels(runs_list)
    plt.xticks(xticks)
    for tick in ax2.xaxis.get_major_ticks():
         tick.label.set_fontsize(10)

    plt.tight_layout()

    margin = (global_space)*len(runs_list)

    plt.legend(fontsize=10)
    #, bbox_to_anchor=(margin, 0), loc='upper left')
    plt.savefig('%s.pdf' % outfile, transparency=True)


def parse_log(tid, runid, fname, tools, tools_id, runs_list):

    nos = -1
    s   = -1

    if (tid, runid) not in tools:
        tools[(tid, runid)] = [{},{}]

    with open(fname) as fp:
        for out in fp.readlines():
            tab = out.strip().split(',')
            
            t = tab[0].split('_')[:2]
            t = '_'.join(t)
            ids = tuple(t.split())
            nos = int(tab[1].split('=')[1])
            s = int(tab[2].split('=')[1])
            if (ids[1]) not in tools[(tid, runid)][0]:
                tools[(tid, runid)][0][ids[1]] = [nos]
                tools[(tid, runid)][1][ids[1]] = [s]
            else:
                tools[(tid, runid)][0][ids[1]].append(nos)
                tools[(tid, runid)][1][ids[1]].append(s)
            runs_list.add(ids[1])
        tools_id.append((tid, runid))
    print (tools_id)


import re
def sorted_nicely( l ):
    """ Sorts the given iterable in the way that is expected.
        Required arguments:
       l -- The iterable to be sorted.
    """
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key = alphanum_key)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-t', '--tool', action='append', dest='tools', default=[], nargs='+')
                        #metavar=('tool_name', 'tool_id', 'nosignal', 'control'))
    parser.add_argument('-o', '--output', action='store', dest='outfile', default='./nosignal',
                        help='Output file where the output files will be stored')
    parser.add_argument('--labels', help='The labels for the plot lines of the ratios.')
    args = parser.parse_args()

    tools_ids = []
    tools = {}
    runs_list = set()
    for xx in args.tools:
        if len(xx) < 3:
            print ("Incorrect number of parameters in curve.")
        if xx[0] not in all_stats:
            print ('ERROR tool %s is not available' % xx[0])
            exit(-1)

        module_ = __import__('tools.' + xx[0].lower(), fromlist=xx[0].title())
        class_ = getattr(module_, xx[0].title())
        operator[xx[0]] = class_()
        parse_log(xx[0], xx[1], xx[2], tools, tools_ids, runs_list)
#        tools[(xx[0], xx[1])] = parse_log(xx[2])
#        tools_ids.append((xx[0], xx[1]))
    print(tools)
    if not os.path.exists(os.path.dirname(args.outfile)):
        try:
            os.makedirs(os.path.dirname(args.outfile))
        except OSError as exc:  # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise

    main(tools_ids, tools, args.labels, sorted_nicely(list(runs_list)), outfile=args.outfile)
