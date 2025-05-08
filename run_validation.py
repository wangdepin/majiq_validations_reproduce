import argparse
import gffutils
import os
import pickle
import numpy as np
from tools import Tools, operator, all_stats
import sys
import matplotlib.pyplot as plt

hash_map = [{}, {}, {}, {}, {}]

__author__    = "Jorge D. Vaquero"
__copyright__ = "Licensed under GNU GPLv3"


class Variation(object):

    def __init__(self, gene_id, transcript, start, end, ratio):
        self.gene = gene_id
        self.transcript_list = [transcript]
        self.ratio = np.copy(ratio)
        self.start = start
        self.end = end

    def update(self, ratio):
        # print(self.ratio, ratio)
        self.ratio += ratio

    def add_transcript(self, transcript):
        self.transcript_list.append(transcript)


def _new_variation(gene_id, transcript, start, end, isexon, ratio, ir=False):
    coord_key = '%s-%s' % (start, end)
    prim_key = '%s:%s' % (gene_id, coord_key)

    if ir:
        prim_key = 'IR:' + prim_key
    indx = 0 if isexon else 1
    try:
        obj = hash_map[indx][prim_key]
        obj.update(ratio)
        obj.add_transcript(transcript)

    except KeyError:
        obj = Variation(gene_id, transcript, start, end, ratio)
        hash_map[indx][prim_key] = obj
        if indx:
            if gene_id in hash_map[2]:
                hash_map[2][gene_id][coord_key] = obj
            else:
                hash_map[2][gene_id] = {coord_key: obj}
    return obj


overlap_g = {}

def _parse_annotation_db(fn, dbfld, ratios):
    db = gffutils.create_db(fn, dbfn='%s/db.tmp' % dbfld, force=True)
    for ii in db.features_of_type('gene'):
        gnid = ii.attributes['ID'][0]
        for txt in db.children(ii, level=1):
            txtid = txt.attributes['ID'][0]
            jstart = -1
            try:
                rat = ratios[txtid]
            except KeyError:
                rat = np.array([0., 0.])
            for ex in db.children(txt, level=1, featuretype='exon', order_by='start'):
                #print gnid, txtid, ex.start, ex.end

                _new_variation(gnid, txtid, ex.start, ex.end, isexon=True, ratio=rat)

                if jstart == -1:
                    jstart = ex.end
                    continue
                jend = ex.start
                _new_variation(gnid, txtid, jstart, jend, isexon=False, ratio=rat)
                jstart = ex.end


def parse_transcript_file(fname, ttable=None, conditions=[]):
    tlb = {}
    gene_tlb = {}
    transcript_length = {}

    cond1 = range(0, len(conditions[0]))
    cond2 = range(len(conditions[0]), len(conditions[1])+len(conditions[0]))

    with open(ttable) as fp:
        for ll in fp.readlines()[1:]:
            if ll.startswith('#'): continue

            tab = ll.strip().split()
            tlb[tab[0]] = tab[1]
            gene_tlb[tab[1]] = tab[2]
            transcript_length[tab[1]] = int(tab[3])
            if tab[2] not in hash_map[3]:
                hash_map[3][tab[2]] = [tab[1]]
            else:
                hash_map[3][tab[2]].append(tab[1])

    ncol = -1

    transcript_count = {}
    gene_total = {}
    with open(fname) as fp:
        for ll in fp.readlines():
            if ll.startswith('#ID'):
                tab = ll.strip().split()
                cond1 = [xid for xid, xx in enumerate(tab[1:]) if xx in conditions[0]]
                cond2 = [xid for xid, xx in enumerate(tab[1:]) if xx in conditions[1]]
                continue

            if ll.startswith('#'):
                continue

            tab = ll.strip().split()
            assert ncol == -1 or ncol == len(tab[1:]), "WRONG transcript expression format"
            ncol = len(tab[1:])
            try:
                trans = tab[0] if ttable is None else tlb[tab[0]]
            except KeyError:
                print ('Error missing transcript: ', tab[0])
                continue
            vals = np.array([float(xx) for xx in tab[1:]])
            values = vals/transcript_length[trans]

            transcript_count[trans] = np.array([np.nanmean(values[cond1]), np.nanmean(values[cond2])])
            where_are_NaNs = np.isnan(transcript_count[trans])
            transcript_count[trans][where_are_NaNs] = 0
            gn = gene_tlb[trans]
            try:
                gene_total[gn] += values
            except KeyError:
                gene_total[gn] = values #np.zeros(shape=ncol)

    return transcript_count


def main_gen_truth(args):

    ratios_dict = parse_transcript_file(args.transcript_expression, args.ttable, [args.cond1, args.cond2])
    tmpdir = '%s/tmp' % args.outDir

    if not os.path.exists(tmpdir):
        os.makedirs(tmpdir)

    _parse_annotation_db(fn=args.annotationDB, dbfld=tmpdir, ratios=ratios_dict)
    hash_map[4] = ratios_dict
    if args.ir_file is not None:
        with open(args.ir_file) as fp:
            for l in fp.readlines():
                if l.startswith('#'): continue
                tab = l.strip().split()
                st, end = tab[1].split('-')
                found = False
                for xx in tab[2].split(';'):
                    try:
                        found = True
                        rat = ratios_dict[xx]
                    except KeyError:
                        print ("WARNING MISSING TRANSCRIPT: ", xx)
                        continue

                    if found:
                        _new_variation(tab[0], xx, st, end, isexon=False, ratio=rat, ir=True)

    with open("%s/annot.pkl" % args.outDir, 'w+b') as fp:
        pickle.dump(hash_map, fp)
    print("Truth File generated in: %s/annot.pkl" % args.outDir)


##

""" VALIDATION PIPELINE"""


def main_validate(args):

    tools = {}
    print(" ".join(sys.argv))

    if args.tool is not None:
        for xx in args.tool:
            if xx[0] not in all_stats:
                print('ERROR tool %s is not available' % xx[0])
                exit(-1)

            module_ = __import__('tools.' + xx[0].lower(), fromlist=xx[0].title())
            class_ = getattr(module_, xx[0].title())
            operator[xx[0]] = class_()
            tools[(xx[0], xx[1])] = xx[2]

    if not os.path.exists(args.outDir):
        os.makedirs(args.outDir)

    with open(args.truth_file, 'rb') as fp:
        hash_map = pickle.load(fp)

    glist = list(hash_map[2].keys())
    if args.gene_list:
        with open(args.gene_list) as fp:
            glist_file = [xx.strip() for xx in fp.readlines()]
        gene_list = [xx for xx in glist if xx in glist_file]
    else:
        gene_list = glist

    stats_dict = {}
    for (tool_id, exec_id), files in tools.items():
        print (tool_id, exec_id, args.ir)
        if args.per_gene:
            stats_dict[(tool_id, exec_id)] = operator[tool_id].validate_chg_genes(files, hash_map,
                                                                                  gene_list,
                                                                                  dpsi_threshold=args.dpsi_t,
                                                                                  pval_threshold=args.pval_t,
                                                                                  ir=args.ir, complex=args.complex,
                                                                                  min_exp=args.min_exp, 
                                                                                  use_score=args.use_score)
        else:
            stats_dict[(tool_id, exec_id)] = operator[tool_id].validate(files, hash_map,
                                                                        dpsi_threshold=args.dpsi_t,
                                                                        pval_threshold=args.pval_t,
                                                                        ir=args.ir, complex=args.complex,
                                                                        min_exp=args.min_exp,
                                                                        use_score=args.use_score)

    with open('%s/dataset_stats.tsv' % args.outDir, 'w+') as fp:
        res_str = Tools.print_header()
        fp.write('%s\n' % res_str)
        print(res_str)
        for tid, dd in stats_dict.items():
            res_str = Tools.print_stats(dd)
            fp.write('%s\t%s\n' % (tid[1], res_str))

            print(tid[1], res_str)

def main_dpsi_delta(args):

    tools = {}
    print(" ".join(sys.argv))

    if args.tool is not None:
        for xx in args.tool:
            if xx[0] not in all_stats:
                print('ERROR tool %s is not available' % xx[0])
                exit(-1)

            module_ = __import__('tools.' + xx[0].lower(), fromlist=xx[0].title())
            class_ = getattr(module_, xx[0].title())
            operator[xx[0]] = class_()
            tools[(xx[0], xx[1])] = xx[2]

    if not os.path.exists(args.outDir):
        os.makedirs(args.outDir)

    with open(args.truth_file, 'rb') as fp:
        hash_map = pickle.load(fp)

    stats_dict = {}
    print("ARGS: " , args)
    for (tool_id, exec_id), files in tools.items():
        stats_dict[(tool_id, exec_id)] = operator[tool_id].dpsi_delta(files, hash_map,
                                                                    dpsi_threshold=args.dpsi_t,
                                                                    ir=False, complex=args.complex,
                                                                    min_exp=args.min_exp)
    plt.figure(figsize=[12, 12])
    plt.rcParams.update({'font.size': 18})
    for xx, yy in stats_dict.items():
        color, lst = operator[xx[0]].get_color(False)
        plt.plot(yy[0], yy[1], label="%s (N=%s)" %(xx[1], yy[2]), c=color, ls=lst)
    plt.legend()
    plt.title(os.path.basename(args.outDir))
    plt.savefig('%s.pdf' % args.outDir)



if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    common = argparse.ArgumentParser(add_help=False)
    common.add_argument('-o', '--output', action='store', dest='outDir', default='./',
                        help='Output directory where the output files will be stored')

    gen_anot = argparse.ArgumentParser(add_help=False)

    gen_anot.add_argument('transcript_expression', type=str,
                          help="True Transcript expression, expected format is as tsv:"
                               "TRANSCRIPT ID\tCOVERAGE[0..Inf)\t..\tCOVERAGE[0..Inf)")
    gen_anot.add_argument('annotationDB', type=str,
                          help="Annotation DB (for example ensembl or refseq) used for the AS tools to quantify RNASeq,"
                               "The transcript IDs should match the transcript_expression file. In case they don't"
                               " match please provide a translation table using --ttable option")
    gen_anot.add_argument('--cond1', nargs='+')
    gen_anot.add_argument('--cond2', nargs='+')
    gen_anot.add_argument('--ttable', action="store", required=True,
                          help="File with the transcript translation between transcript_Expression file and annotation "
                               "DB")
    gen_anot.add_argument('--ir-file', action="store", help="This file of intron retention")

    validate = argparse.ArgumentParser(add_help=False)
    validate.add_argument('truth_file', type=str, help='Pickle file created by gen_truth pipeline.')
    validate.add_argument('-t', '--tool', action='append', nargs=3, metavar=('tool_name', 'exec_id', 'file1'))
    validate.add_argument('-d', '--dpsi_threshold', action='store', dest='dpsi_t', default=0.2, type=float,
                          help='Threshold for minimum expected DPSI to be considered true changing')
    validate.add_argument('-p', '--pval_threshold', action='store', dest='pval_t', default=0.05, type=float,
                          help='Threshold for maximum pval to be considered true changing')
    validate.add_argument('--permissive', action="store_true", default=False,
                          help='Adding permivissive flag will disable accepting significant cases (pval<=pval_thesh) '
                               'with low deltapsi as FN')
    validate.add_argument('--ir', action="store_true", default=False)
    validate.add_argument('--per-gene', action="store_true", default=False)
    validate.add_argument('--use-score', action='store_true', default=False)
    validate.add_argument('--complex', action='store_true', default=False,
                          help='Use complex event in whippet')
    validate.add_argument('--min-exp', action='store', dest='min_exp', default=0.5, type=float,
                          help='Min experiments with CI>= 0.10  for whippet')
    validate.add_argument('--gene-list', action="store", help="This option defines a file with a list of genes, if "
                                                              "this is defined, we will validate only those genes. This"
                                                              " option has no effect on per_event validation")

    dpsi_delta = argparse.ArgumentParser(add_help=False)
    dpsi_delta.add_argument('truth_file', type=str, help='Pickle file created by gen_truth pipeline.')
    dpsi_delta.add_argument('-t', '--tool', action='append', nargs=3, metavar=('tool_name', 'exec_id', 'file1'))
    dpsi_delta.add_argument('-d', '--dpsi_threshold', action='store', dest='dpsi_t', default=-1, type=float,
                          help='Threshold for minimum expected DPSI to be considered true changing')
    dpsi_delta.add_argument('--complex', action='store_true', default=False,
                        help='Use complex event in whippet')
    dpsi_delta.add_argument('--min-exp', action='store', dest='min_exp', default=0.5, type=float,
                        help='Min experiments with CI>= 0.10 for whippet')


    subparsers = parser.add_subparsers(help='')
    parser_preprocess = subparsers.add_parser('gen_truth',
                                              help='Generates the truth based on transcript expression '
                                                   'values and the annotationDB used to base the '
                                                   'transcript expression file',
                                              parents=[common, gen_anot])
    parser_preprocess.set_defaults(func=main_gen_truth)

    parser_validate = subparsers.add_parser('validate', help='Validates the tool outputs using the truth generated '
                                                             'by gen_truth',
                                            parents=[common, validate])
    parser_validate.set_defaults(func=main_validate)
    parser_dpsi_delta = subparsers.add_parser('dpsi-delta', help='Generates a CDF of the difference between '
                                                                 'deltapsi and groundtruth dpsi',
                                            parents=[common, dpsi_delta])
    parser_dpsi_delta.set_defaults(func=main_dpsi_delta)
    args = parser.parse_args()

    args.func(args)

