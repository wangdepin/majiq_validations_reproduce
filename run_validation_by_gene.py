import argparse
import gffutils
import os
import pickle
import numpy as np
from tools import Tools, operator, all_stats
import sys

hash_map = [{}, {}]

__author__    = "Jorge D. Vaquero"
__copyright__ = "Licensed under GNU GPLv3"

#
# class Variation(object):
#
#     def __init__(self, gene_id, transcript, start, end, ratio):
#         self.gene = gene_id
#         self.transcript_list = [transcript]
#         self.ratio = np.copy(ratio)
#         self.start = start
#         self.end = end
#
#     def update(self, ratio):
#         # print(self.ratio, ratio)
#         self.ratio += ratio
#
#     def add_transcript(self, transcript):
#         self.transcript_list.append(transcript)

#
# def _new_variation(gene_id, transcript, start, end, isexon, ratio, ir=False):
#     prim_key = '%s:%s-%s' % (gene_id, start, end)
#     if ir:
#         prim_key = 'IR:' + prim_key
#     indx = 0 if isexon else 1
#     try:
#         obj = hash_map[indx][prim_key]
#         obj.update(ratio)
#         obj.add_transcript(transcript)
#     except KeyError:
#         obj = Variation(gene_id, transcript, start, end, ratio)
#         hash_map[indx][prim_key] = obj
#
#     return obj
#
#
# overlap_g = {}
#
# def _parse_annotation_db(fn, dbfld, ratios):
#     db = gffutils.create_db(fn, dbfn='%s/db.tmp' % dbfld, force=True)
#     for ii in db.features_of_type('gene'):
#         gnid = ii.attributes['ID'][0]
#         for txt in db.children(ii, level=1):
#             txtid = txt.attributes['ID'][0]
#             jstart = -1
#             try:
#                 rat = ratios[txtid]
#             except KeyError:
#                 rat = np.array([0., 0.])
#             for ex in db.children(txt, level=1, featuretype='exon', order_by='start'):
#                 #print gnid, txtid, ex.start, ex.end
#
#                 _new_variation(gnid, txtid, ex.start, ex.end, isexon=True, ratio=rat)
#
#                 if jstart == -1:
#                     jstart = ex.end
#                     continue
#                 jend = ex.start
#                 _new_variation(gnid, txtid, jstart, jend, isexon=False, ratio=rat)
#                 jstart = ex.end
#

def find_changing_genes_by_isoform(fname, ttable, group1, group2):
    tlb = {}
    gene_tlb = {}
    transcript_length = {}

    cond1 = range(0, len(group1))
    cond2 = range(len(group1), len(group2)+len(group1))

    with open(ttable) as fp:
        for ll in fp.readlines()[1:]:
            if ll.startswith('#') : continue

            tab = ll.strip().split()
            tlb[tab[0]] = tab[1]
            gene_tlb[tab[1]] = tab[2]
            transcript_length[tab[1]] = int(tab[3])

    ncol = -1

    gene_max_change = {}
    gene_total = {}
    transcript_count = {}
    with open(fname) as fp:
        for ll in fp.readlines():
            if ll.startswith('#ID'):
                tab = ll.strip().split()
                cond1 = [xid for xid, xx in enumerate(tab[1:]) if xx in group1]
                cond2 = [xid for xid, xx in enumerate(tab[1:]) if xx in group2]
                continue

            if ll.startswith('#'):
                continue

            tab = ll.strip().split()
            assert ncol == -1 or ncol == len(tab[1:]), "WRONG transcript expression format"
            ncol = len(tab[1:])
            try:
                trans = tab[0] if ttable is None else tlb[tab[0]]
            except KeyError:
                continue
            values = np.array([float(xx) for xx in tab[1:]])
            values /= transcript_length[trans]
            transcript_count[trans] = np.array([np.nanmean(values[cond1]), np.nanmean(values[cond2])])
            gn = gene_tlb[trans]
            try:
                gene_total[gn] += transcript_count[trans]
            except KeyError:
                gene_total[gn] = transcript_count[trans]

    for trans, (cnt1, cnt2) in transcript_count.items():
        gn = gene_tlb[trans]
        psi1 = cnt1 / (gene_total[gn][0])
        psi2 = cnt2 / (gene_total[gn][1])
        dpsi = abs(psi1 - psi2)
#        print(gn, cnt1, cnt2, psi1, psi2,"gene_total", gene_total[gn][0], gene_total[gn][1], dpsi)
        if (np.isnan(dpsi)): 
            continue
        try:
            gene_max_change[gn] = max(gene_max_change[gn], dpsi)
        except KeyError:
            gene_max_change[gn] = dpsi
    return gene_max_change


def main_all(args):

    tools = {}
    if args.tool is not None:
        for xx in args.tool:
            if xx[0] not in all_stats:
                print ('ERROR tool %s is not available' % xx[0])
                exit(-1)

            module_ = __import__('tools.' + xx[0].lower(), fromlist=xx[0].title())
            class_ = getattr(module_, xx[0].title())
            operator[xx[0]] = class_()
            tools[(xx[0], xx[1])] = xx[2]

    if not os.path.exists(args.outDir):
        os.makedirs(args.outDir)

    ratios_dict = find_changing_genes_by_isoform(args.transcript_expression, args.ttable, args.cond1, args.cond2)

    with open("%s/gene_change.pkl" % args.outDir, 'w+b') as fp:
        pickle.dump(ratios_dict, fp)

    print(" ".join(sys.argv))


    stats_dict = {}
    for (tool_id, exec_id), files in tools.items():
        print(tool_id, exec_id, args.ir)
        stats_dict[(tool_id, exec_id)] = operator[tool_id].validate_chg_genes(files, ratios_dict,
                                                                    dpsi_threshold=args.dpsi_t,
                                                                    pval_threshold=args.pval_t,
                                                                    permissive=args.permissive,
                                                                    ir=args.ir, complex=args.complex,
                                                                    min_exp=args.min_exp)

    with open('%s/dataset_stats.tsv' % args.outDir, 'w+') as fp:
        res_str = Tools.print_header()
        fp.write('%s\n' % res_str)
        print (res_str)
        for tid, dd in stats_dict.items():
            res_str = Tools.print_stats(dd)
            fp.write('%s\t%s\n' % (tid[1], res_str))

            print (tid[1], res_str)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--output', action='store', dest='outDir', default='./',
                        help='Output directory where the output files will be stored')

    parser.add_argument('transcript_expression', type=str,
                         help="True Transcript expression, expected format is as tsv:"
                              "TRANSCRIPT ID\tCOVERAGE[0..Inf)\t..\tCOVERAGE[0..Inf)")

    parser.add_argument('--cond1', nargs='+')
    parser.add_argument('--cond2', nargs='+')
    parser.add_argument('--ttable', action="store", required=True,
                          help="File with the transcript translation between transcript_Expression file and annotation "
                               "DB")


    parser.add_argument('-t', '--tool', action='append', nargs=3, metavar=('tool_name', 'exec_id', 'file1'))
    parser.add_argument('-d', '--dpsi_threshold', action='store', dest='dpsi_t', default=0.2, type=float,
                          help='Threshold for minimum expected DPSI to be considered true changing')
    parser.add_argument('-p', '--pval_threshold', action='store', dest='pval_t', default=0.05, type=float,
                          help='Threshold for maximum pval to be considered true changing')
    parser.add_argument('--permissive', action="store_true", default=False,
                          help='Adding permivissive flag will disable accepting significant cases (pval<=pval_thesh) '
                               'with low deltapsi as FN')
    parser.add_argument('--ir', action="store_true", default=False)
    parser.add_argument('--use-score', action='store_true', default=False,
                        help='Use tnom score, only for TNOM in majiq het')

    #whippet
    parser.add_argument('--complex', action='store_true', default=False,
                        help='Use complex event in whippet')
    parser.add_argument('--min-exp', action='store', dest='min_exp', default=1.0, type=float,
                        help='Min experiments with CI>= 0.10')

    args = parser.parse_args()

    main_all(args)

