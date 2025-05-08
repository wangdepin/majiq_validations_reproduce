from os.path import dirname, basename, isfile
import glob
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
from statsmodels.distributions.empirical_distribution import ECDF
import math

__author__    = "Jorge D. Vaquero"
__copyright__ = "Licensed under GNU GPLv3"

modules = glob.glob(dirname(__file__)+"/*.py")
all_stats = [basename(f)[:-3] for f in modules if isfile(f) if basename(f)[:-3] not in ['__init__', 'skeleton']]
operator = {}


class ToolsFactory:

    factories = {}

    @staticmethod
    def add_factory(id_stat, tools_factory):
        ToolsFactory.factories.put[id] = tools_factory

    # A Template Method:
    @staticmethod
    def create_tools(id_tool):
        if not ToolsFactory.factories.has_key(id_tool):
            ToolsFactory.factories[id_tool] = eval(id_tool + '.Factory()')
        return ToolsFactory.factories[id_tool].create()


class Tools(object):

    index_vals = {'dpsi': 1, 'pval': 2}
    lstyles = ['solid', 'dashed', 'dotted', 'dashdot']
    idx_stats = {'NUM_EVENTS': 0, 'TOTAL_FOUND': 1, 'NOT_FOUND': 2,  'GREY_EVENTS': 3, 'TOTAL_CHG': 4, 'TOTAL_NO_CHG': 5,
                 'TP': 6, 'TN': 7, 'FP': 8, 'FN': 9}
    key_stats = ['NUM_EVENTS', 'TOTAL_FOUND', 'NOT_FOUND', 'GREY_EVENTS', 'TOTAL_CHG', 'TOTAL_NO_CHG',
                 'TP', 'TN', 'FP', 'FN']

    @staticmethod
    def print_stats(res_dic):
        rs = ''
        for kk in Tools.key_stats:
            rs += '%s\t' % res_dic[Tools.idx_stats[kk]]
        try:
    
            tp =  float(res_dic[Tools.idx_stats['TP']])
            fp = float(res_dic[Tools.idx_stats['FP']])
            tn = float(res_dic[Tools.idx_stats['TN']])
            fn = float(res_dic[Tools.idx_stats['FN']])
    
    
            sens = tp / (tp + fn)
            spec = tn / (tn + fp)
            prec = tp / (tp + fp)
            fdr =  fp / (tp + fp)
            fnr = fn / (tp + fn)
            # fpr = float(res_dic[Tools.idx_stats['FP']]) / (res_dic[Tools.idx_stats['TN']] + res_dic[Tools.idx_stats['FP']])
            fscore = (2*prec*sens) / (prec+sens)
            mcc = ((tp*tn) - (fp*fn)) / math.sqrt((tn+fn)*(fp+tp)*(tn+fp)*(fn+tp))
    
        except ZeroDivisionError:
            print("TP: ", res_dic[Tools.idx_stats['TP']], " FN: ", res_dic[Tools.idx_stats['FN']],
                  " TN: ", res_dic[Tools.idx_stats['TN']], " FP: ", res_dic[Tools.idx_stats['FP']])
            sens = 'NA'
            spec = 'NA'
            fdr  = 'NA'
            fnr  = 'NA'
            fscore = 'NA'
            mcc  = 'NA'

            #raise

        rs += '%s\t%s\t%s\t%s\t%s\t%s' % (sens, spec, fdr, fnr, fscore, mcc)
        return rs

    @staticmethod
    def print_header():
        rs = 'Tool\t'
        for kk in Tools.key_stats:
            rs += '%s\t' % kk
        rs += "Sensitivity\tSpecificity\tFDR\tFNR\tFSCORE\tMCC"

        return rs

    @staticmethod
    def _psi(a, b):
        try:
            c = np.divide(a, b)
            # print("IN PSI,", a, b, c)
            return abs(c)
        except ZeroDivisionError:
            return 0.

    @staticmethod
    def _dpsi(a, b):
        try:
            c = np.divide(a, b)
            return abs(c[0] - c[1])
        except ZeroDivisionError:
            return -1
    @staticmethod
    def _sort_rank(rank, pval=False):
        if pval:
            print("Sorting by pval")
            #return sorted(rank, key=lambda x: (-x[Tools.index_vals['pval']], abs(x[Tools.index_vals['dpsi']])),
            return sorted(rank, key=lambda x: (-x[Tools.index_vals['pval']]), reverse=True)
        else:
            return sorted(rank, key=lambda x: (abs(x[Tools.index_vals['dpsi']]), -x[Tools.index_vals['pval']]),
                                reverse=True)

    @staticmethod
    def _print_stats(nexp, nfdr, nexp_fdr, dpsi_thresh, pval_thresh, method):
        print ("[%s] #FDR < %.2f: %d" % (method, pval_thresh, nfdr))
        print ("[%s] #E(Delta(PSI))>%.2f: %d" % (method, dpsi_thresh,  nexp))
        print ("[%s] #E(Delta(PSI))>%.2f and FDR<%.2f: %d" % (method, dpsi_thresh, pval_thresh, nexp_fdr))

    @staticmethod
    def _rr_vals(rank1, rank2, method=""):
        # max_events = len([xx for xx in rank1 if xx[-1]])
        max_events = len(rank1)
        LSVs_1 = list(map(lambda xx: xx[0], rank1))
        LSVs_2 = list(map(lambda xx: xx[0], rank2))
        #    frac_max_N = 1.0/max_N
        ratios = np.array([np.intersect1d(LSVs_1[:n + 1], LSVs_2[:n + 1]).size * (1.0 / (n + 1))
                           for n in range(max_events)])
        print("[%s] N=%s, RR=%s" % (method, max_events, ratios[max_events - 1]))
        return ratios, max_events

    @classmethod
    def validate_chg_genes(method, outputfile, validation_dict, gene_list, dpsi_threshold=0.2, pval_threshold=0.05, **kargs):
        values = [0] * len(method.idx_stats)
        gene_dict = method._validate_chg_genes(outputfile, validation_dict, values, **kargs)
        print(len(gene_list))
        for gn in gene_list:
            values[method.idx_stats['NUM_EVENTS']] += 1
            if gn in gene_dict:
                values[method.idx_stats['TOTAL_FOUND']] += 1
            else:
                values[method.idx_stats['NOT_FOUND']] += 1
                continue

            tp = 0
            fp = 0 
            tn = 0 
            fn = 0
            gn_changing = False

            for ev_vals in gene_dict[gn]:
                ratio, dpsi, probs = ev_vals 
                if ratio >= dpsi_threshold:
                    '''True Positive'''
                    gn_changing |= True
                    if abs(dpsi) >= dpsi_threshold and probs <= pval_threshold :
                        tp += 1

                # else:
                    elif abs(dpsi) < 0.05 and probs > pval_threshold:
                        fn += 1
                        print("FN", gn, "ratio =", ratio, " dpsi = ", dpsi, '>=',
                              dpsi_threshold, ':: Prob=', probs, '>=', pval_threshold)

                elif ratio < 0.05:
                    '''True Negative'''
                    gn_changing |= False

                    if abs(dpsi) >= dpsi_threshold and probs <= pval_threshold:
                        fp += 1
                        print("FP", gn, " ratio = ", ratio, " dpsi = ", abs(dpsi), '>=',
                              dpsi_threshold, ':: Prob=', probs, '>=', pval_threshold)

                    elif abs(dpsi) < 0.05 and probs > pval_threshold:
                        tn += 1

            if gn_changing:
                values[method.idx_stats['TOTAL_CHG']] += 1
            else:
                values[method.idx_stats['TOTAL_NO_CHG']] += 1

            if fp > 0:
                values[method.idx_stats['FP']] += 1
            elif fn > 0:
                values[method.idx_stats['FN']] += 1
            elif tp > 0:
                values[method.idx_stats['TP']] += 1
            elif tn > 0:
                values[method.idx_stats['TN']] += 1
                    
        return values
    
    @classmethod
    def validate(method, outputfile, validation_dict, dpsi_threshold=0.2, pval_threshold=0.05, **kargs):
        values = [0] * len(method.idx_stats)
        events = method._validate(outputfile, validation_dict, values, **kargs)
        for ev, ev_vals in events.items():
            values[method.idx_stats['NUM_EVENTS']] += 1
            ratio, dpsi, probs = ev_vals 
            if ratio >= dpsi_threshold:
                '''True Positive'''
                values[method.idx_stats['TOTAL_CHG']] += 1

                if abs(dpsi) >= dpsi_threshold and probs <= pval_threshold :
                    idx = method.idx_stats['TP']
                    values[idx] += 1

                # else:
                elif abs(dpsi) < 0.05 and probs > pval_threshold:
                    print("FN", ev, "ratio =", ratio, " dpsi = ", dpsi, '>=',
                          dpsi_threshold, ':: Prob=', probs, '>=', pval_threshold)
                    idx = method.idx_stats['FN']
                    values[idx] += 1

            elif ratio <= 0.05:
                '''True Negative'''
                values[method.idx_stats['TOTAL_NO_CHG']] += 1

                if abs(dpsi) >= dpsi_threshold and probs <= pval_threshold:
                    print("FP", ev, " ratio = ", ratio, " dpsi = ", abs(dpsi), '>=',
                          dpsi_threshold, ':: Prob=', probs, '>=', pval_threshold)
                    idx = method.idx_stats['FP']
                    values[idx] += 1

                elif abs(dpsi) < 0.05 and probs > pval_threshold:
                    idx = method.idx_stats['TN']
                    values[idx] += 1

            else:
                values[method.idx_stats['GREY_EVENTS']] += 1
                continue
        return values

    @classmethod
    def dpsi_delta(method, files, validation_dict, dpsi_threshold=-1, **kargs):
        values = method._dpsi_delta(files, validation_dict, dpsi_threshold, **kargs)
        ecdf = ECDF(values)
#        print(" @ ", values, len(values))
        return ecdf.x, ecdf.y, len(values)

    @classmethod
    def iir_calc(method, file_prefix1, file_prefix2, dpsi_thresh=0.2, pval_thresh=0.05, msg='', **kargs):
        n_nosignal = method._count_sig_events(file_prefix1, dpsi_thresh=dpsi_thresh, pval_thresh=pval_thresh, **kargs)
        n_signal = method._count_sig_events(file_prefix2, dpsi_thresh=dpsi_thresh, pval_thresh=pval_thresh, **kargs)
        if n_signal > 0:
            iir = float(n_nosignal) / n_signal
        else:
            n_signal = -1; n_nosignal = -1; iir = -0.1
        print ("%s, N_nosignal=%s, N_signal=%s, IIR=%s" % (msg, n_nosignal, n_signal, iir))

    @staticmethod
    def dump_values_corr(data1, data2, output, labels=None):
        with open("%s.tsv" % output, "w+") as ofp:
            for i in range(len(data1)):
                if labels is not None:
                    s = "%s\t" % labels[i]
                else:
                    s = ""
                s += "%s\t%s" % (data1[i], data2[i])
                ofp.write("%s\n" % s);

    @staticmethod
    def plot_corr(data1, data2, title='Title', xax='RTPCR', yax='Tool', f_idx=0, name='rtpcr',
                  lims=(-1, 1), cline='#262626'):
        plt.figure(f_idx)
        fig, ax = plt.subplots(1)
        ax.set_aspect('equal')
        plt.scatter(data1, data2,
                    facecolors=cline, edgecolors=cline,
                    alpha=0.95)
        plt.axvline(0, c='#262626', linewidth=1.5, alpha=0.9)
        plt.axhline(0, c='#262626', linewidth=1.5, alpha=0.9)
        plt.ylim(lims)
        plt.xlim(lims)
        fit = Tools.plot_fit(data1, data2)
        plt.plot(data1, fit(data1), '#262626', linewidth=1.5, alpha=0.9)
        plt.xlabel(xax)
        plt.ylabel(yax)
        plt.title(title)
        plt.savefig('%s.pdf' % name, box_inches='tight')

    @staticmethod
    def plot_fit(x_data, y_data, return_stats=False):
        r = scipy.stats.pearsonr(x_data, y_data)
        fit = np.polyfit(x_data, y_data, 1)
        fit_fn = np.poly1d(fit)
        plt.plot(x_data, fit_fn(x_data), '#262626', linewidth=1.0, alpha=0.7)
        print ("R2: %s" % (r[0] * r[0]))
        plt.annotate('Pearson\' r:\n%.3f\nR^2:%.3f' % (r[0], r[0] * r[0]), fontsize=10,
                   xy=(0.15, 0.9), xycoords='axes fraction')
        if return_stats == False:
            print (r, r[0] * r[0])
            return fit_fn
        if return_stats == True:
            return r[0], r[0] * r[0], r[1]


