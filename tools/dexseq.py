from tools import Tools
import brewer2mpl
import numpy as np

__author__    = "Jorge D. Vaquero"
__copyright__ = "Licensed under GNU GPLv3"

class Dexseq(Tools):
    __pval_cache__ = {}

    plot_count = 0
    colors = np.array(brewer2mpl.get_map('Paired', 'qualitative', 11).mpl_colors)[[7, 6]]
    lstyle_list = ['solid', 'dashed', 'dotted']

    class Factory:
        def create(self): return Dexseq()

    """private functions"""
    @staticmethod
    def _rank_dexseq(dexseq_file, dpsi_thresh=0.2, pval_thresh=0.05,
                     pval_f=False, step1_list=None):
        rank = []
        dexseq_nn = 0
        list_events = []
        for line in open(dexseq_file).readlines()[1:]:

            sline = line.strip().split()
            if 'NA' in sline:
                    continue
            log2fold = float(sline[9])
            padj = float(sline[6])
            ev_name = sline[0]
            pass_thres = int(padj < pval_thresh and log2fold > 4)
            if step1_list is None:
                list_events.append(ev_name)
            elif ev_name not in step1_list:
                continue
            dexseq_nn += pass_thres
            if padj <= pval_thresh:
                rank.append([ev_name, log2fold, padj, padj <= pval_thresh])

        Dexseq._print_stats(rank, dpsi_thresh=dpsi_thresh, pval_thresh=pval_thresh, method='DexSeq')
        rank = Dexseq._sort_rank(rank, pval=pval_f)
        return rank, list_events

    """Public functions"""
    @staticmethod
    def rr_rank(file1, file2, dpsi_thresh=0.2, pval_thresh=0.05, pval=False, use_score=False):
        rank1, ev_list = Dexseq._rank_dexseq(file1, dpsi_thresh=dpsi_thresh,
                                             pval_thresh=pval_thresh, pval_f=pval)
        rank2, nev2 = Dexseq._rank_dexseq(file2, dpsi_thresh=dpsi_thresh,
                                          pval_thresh=pval_thresh,
                                          step1_list=ev_list, pval_f=pval)

        ratios = Dexseq._rr_vals(rank1, rank2, method='DexSeq')

        return ratios

    @staticmethod
    def get_color(subsample=False):

        if subsample:
            color = Dexseq.colors[1]
            lstyle = Dexseq.lstyle_list[Dexseq.plot_count % len(Dexseq.lstyle_list)]
        else:
            color = Dexseq.colors[Dexseq.plot_count % len(Dexseq.colors)]
            lstyle = '-'

        Dexseq.plot_count += 1

        return color, lstyle