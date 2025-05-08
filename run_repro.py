import argparse
from tools import operator, all_stats
import pickle
import os

__author__    = "Jorge D. Vaquero"
__copyright__ = "Licensed under GNU GPLv3"


def main(tools, pval=0.05, dpsi=0.2, outDir="./", sort_pval=False, ir=False, only_pval=False, use_score=False,
             complex=False, min_exp=1.0, pval_median=False, use_overlap=False):

    for method, files in tools.items():
        ratios = operator[method].rr_rank(files[0], files[1], dpsi_thresh=dpsi,
                                          pval_thresh=pval, pval=sort_pval, ir=ir, only_pval=only_pval,
                                          use_score=use_score, complex=complex, min_exp=min_exp, pval_median=pval_median, use_overlap=use_overlap)
        with open("%s/%s.ratios.pickle" % (outDir, method), 'w+b') as fp:
            pickle.dump(ratios, fp)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-t', '--tool', action='append', nargs=3, metavar=('tool_name', 'file1', 'file2'))
    parser.add_argument('-o', '--output', action='store', dest='outDir', default='./',
                        help='Output directory where the output files will be stored')
    parser.add_argument('-d', '--dpsi_threshold', action='store', dest='dpsi_t', default=0.2, type=float,
                        help='Threshold for minimum expected DPSI to be considered true changing')
    parser.add_argument('-p', '--pval_threshold', action='store', dest='pval_t', default=0.05, type=float,
                        help='Threshold for maximum pval to be considered true changing')
    parser.add_argument('--s_pval', action='store_true', default=False,
                        help='Sort by pval and dpsi')
    parser.add_argument('--only_pval', action='store_true', default=False,
                        help='filter only by pval')
    parser.add_argument('--enable_ir', action='store_true', default=False,
                        help='Use IR for reproducibility')
    parser.add_argument('--use-score', action='store_true', default=False,
                        help='Use tnom score, only for TNOM in majiq het')
    parser.add_argument('--use-overlap', action='store_true', default=False,
                        help='Use overlaping lsv for RR, default: False')
    parser.add_argument('--pval-median', action='store_true', default=False,
                        help='Use pvalue median instead of pval sampling in majiq het')

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
        tools[xx[0]] = xx[1:]

    if not os.path.exists(args.outDir):
        os.makedirs(args.outDir)

    main(tools, pval=args.pval_t, dpsi=args.dpsi_t, outDir=args.outDir, sort_pval=args.s_pval, ir=args.enable_ir,
         only_pval=args.only_pval, use_score=args.use_score, complex=args.complex, min_exp=args.min_exp,
         pval_median=args.pval_median, use_overlap=args.use_overlap)
