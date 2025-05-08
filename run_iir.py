import argparse
from tools import operator, all_stats
import pickle
import os

__author__    = "Jorge D. Vaquero"
__copyright__ = "Licensed under GNU GPLv3"


def main(tools, pval=0.05, dpsi=0.2, runid='', onlypval=False, use_score=False, use_overlap=False):

    for method, files in tools.items():
        operator[method[0]].iir_calc(files[0], files[1], dpsi_thresh=dpsi, pval_thresh=pval,
                                     msg="%s %s" %(method[0], method[1]), onlypval=onlypval, use_score=use_score,
                                     complex=args.complex, min_exp=args.min_exp, use_overlap=use_overlap)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-t', '--tool', action='append', nargs=4, metavar=('tool_name', 'comp_id', 'file1', 'file2'))
    parser.add_argument('-m', '--run_id', action='store', dest='runid', default='',
                        help='Output directory where the output files will be stored')
    parser.add_argument('-d', '--dpsi_threshold', action='store', dest='dpsi_t', default=0.2, type=float,
                        help='Threshold for minimum expected DPSI to be considered true changing')
    parser.add_argument('-p', '--pval_threshold', action='store', dest='pval_t', default=0.05, type=float,
                        help='Threshold for maximum pval to be considered true changing')
    parser.add_argument('--only_pval', action='store_true', default=False,
                        help='filter only by pval')
    parser.add_argument('--use-score', action='store_true', default=False,
                        help='Use tnom score, only for TNOM in majiq het')
    parser.add_argument('--use-overlap', action='store_true', default=False,
                        help='Use overlaping lsv for RR, default: False')
    #whippet
    parser.add_argument('--complex', action='store_true', default=False,
                        help='Use complex event in whippet')
    parser.add_argument('--min-exp', action='store', dest='min_exp', default=0.5, type=float,
                        help='Min experiments with CI>= 0.10')

    args = parser.parse_args()

    tools = {}

    for xx in args.tool:
        if xx[0] not in all_stats:
            print ('ERROR tool %s is not available' % xx[0])
            exit(-1)

        module_ = __import__('tools.' + xx[0].lower(), fromlist=xx[0].title())
        class_ = getattr(module_, xx[0].title())
        operator[xx[0]] = class_()
        tools[(xx[0], xx[1])] = xx[2:]

    main(tools, pval=args.pval_t, dpsi=args.dpsi_t, runid=args.runid, onlypval=args.only_pval, use_score=args.use_score, use_overlap=args.use_overlap)
