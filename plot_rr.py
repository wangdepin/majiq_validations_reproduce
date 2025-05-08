import matplotlib
matplotlib.use('agg')
import argparse
from tools import operator, all_stats
import pickle
import os
import numpy as np
import matplotlib.pyplot as plt
import errno

__author__    = "Jorge D. Vaquero"
__copyright__ = "Licensed under GNU GPLv3"

lstyles = ['solid', 'dashed', 'dotted', 'dashdot']

def main_plot_rr(curves_ids, curves, outfile='./rr_plot.pdf', max_N=2000, subsamples=False):

    plt.figure(figsize=[16, 16])

    for (tool, tid) in curves_ids:
        files = curves[(tool, tid)]
        ratios = np.zeros(shape=(len(files), max_N), dtype=float)
        cap_n = []
        try:
            for fidx, ff in enumerate(files):
                with open(ff, 'rb') as fp:
                    vals = pickle.load(fp)[0]
                    n = min(vals.shape[0], max_N)
                    ratios[fidx, :n] = vals[:n]
                    cap_n.append(n)

            ratios[ratios == 0] = np.nan
            mean_ratios = np.nanmean(ratios, axis=0)
            mu_n = int(np.mean(cap_n))
            std_ratios = np.nanstd(ratios, axis=0)
            color, lst = operator[tool].get_color(subsamples)
            plot_lab = '%s (N = %d,\nRR=%.2f)' % (tid, mu_n, mean_ratios[mu_n-1])

            print (mu_n, mean_ratios[:mu_n], std_ratios[:mu_n])
            plt.plot(np.arange(mu_n), mean_ratios[:mu_n], lw=6, c=color, label=plot_lab, ls=lst)
            plt.fill_between(np.arange(mu_n), mean_ratios[:mu_n] - std_ratios[:mu_n],
                             mean_ratios[:mu_n] + std_ratios[:mu_n],
                             facecolor=color, alpha=.25)
        except FileNotFoundError as e:
            print(e)
            continue

    name = os.path.basename(outfile)
    plt.xlabel('Number of events detected')
    plt.ylabel('Fraction of events reproduced')
    plt.title(name)
    plt.legend(loc=0)
    plt.xlim(0, max_N)
    plt.ylim(0, 1)
    plt.tight_layout()
    plt.savefig('%s.pdf' % outfile, transparency = True)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-c', '--curve', action='append', dest='curve', default=[], nargs='+')#, metavar=('tool_name', 'id', 'list_of_files'))
    parser.add_argument('-o', '--output', action='store', dest='outfile', default='./',
                        help='Output file where the output files will be stored')
    parser.add_argument('-m', '--max_N', action='store', dest='max_n', default=2000, type=int)
    parser.add_argument('--subsamples', action='store_true', dest='subsamples', default=False)
    args = parser.parse_args()

    curves_d = {}

    curves_ids = []
    for xx in args.curve:
        if len(xx) < 3:
            print ("Incorrect number of parameters in curve.")
        if xx[0] not in all_stats:
            print ('ERROR tool %s is not available' % xx[0])
            exit(-1)

        module_ = __import__('tools.' + xx[0].lower(), fromlist=xx[0].title())
        class_ = getattr(module_, xx[0].title())
        operator[xx[0]] = class_()
        if (xx[0], xx[1]) in curves_d:
            print ("Two curves has the same tool and ID, (tool, ID) should be unambiguos")
        curves_d[(xx[0], xx[1])] = xx[2:]
        curves_ids.append((xx[0], xx[1]))
    dirname = os.path.dirname(args.outfile) or '.'
    if not os.path.exists(dirname):
        try:
            os.makedirs(dirname, exist_ok=True)
        except OSError as exc:  # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise

    main_plot_rr(curves_ids, curves=curves_d, outfile=args.outfile, max_N=args.max_n, subsamples=args.subsamples)
