from tools import Tools
import brewer2mpl
import numpy as np

__author__    = "Jorge D. Vaquero"
__copyright__ = "Licensed under GNU GPLv3"

class Rmats(Tools):
    __pval_cache__ = {}
    col_vals = {'SE':   {'last_coord': 10, 'pval': 18, 'FDR': 19, 'dPSI': 22},
                'A3SS': {'last_coord': 10, 'pval': 18, 'FDR': 19, 'dPSI': 22},
                'A5SS': {'last_coord': 10, 'pval': 18, 'FDR': 19, 'dPSI': 22},
                'MXE':  {'last_coord': 12, 'pval': 20, 'FDR': 21, 'dPSI': 24},
                'RI':   {'last_coord': 10}}

    plot_count = 0
    colors = np.array(brewer2mpl.get_map('Paired', 'qualitative', 11).mpl_colors)[[2, 3]]
    lstyle_list = ['solid', 'dashed', 'dotted']

    class Factory:
        def create(self): return Rmats()

    """private functions"""
    @staticmethod
    def __SE_ratios(validation_dict, coords, gene_id, strand=None):

        ex_start, ex_end = int(coords[0]) + 1, coords[1]
        upstream_start, upstream_end = int(coords[2]) + 1, coords[3]
        downstream_start, downstream_end = int(coords[4]) + 1, coords[5]

        ex_vals = validation_dict['%s:%s-%s' % (gene_id, ex_start, ex_end)]
        upstream = validation_dict['%s:%s-%s' % (gene_id, upstream_start, upstream_end)]
        downstream = validation_dict['%s:%s-%s' % (gene_id, downstream_start, downstream_end)]

        ratio1 = Rmats._dpsi(ex_vals.ratio, upstream.ratio + ex_vals.ratio)
        ratio2 = Rmats._dpsi(ex_vals.ratio, downstream.ratio + ex_vals.ratio)

        return max(ratio1, ratio2)

    @staticmethod
    def __MXE_ratios(validation_dict, coords, gene_id, strand=None):
        ex1_start, ex1_end = int(coords[0]) + 1, coords[1]
        ex2_start, ex2_end = int(coords[2]) + 1, coords[3]

        ex1_vals = validation_dict['%s:%s-%s' % (gene_id, ex1_start, ex1_end)]
        ex2_vals = validation_dict['%s:%s-%s' % (gene_id, ex2_start, ex2_end)]

        ratio = Rmats._dpsi(ex1_vals.ratio, ex1_vals.ratio + ex2_vals.ratio)

        return ratio

    @staticmethod
    def __A3SS_ratios(validation_dict, coords, gene_id, strand=None):
        ex1 = validation_dict['%s:%s-%s' % (gene_id, int(coords[0]) + 1, coords[1])]
        ex2 = validation_dict['%s:%s-%s' % (gene_id, int(coords[2]) + 1, coords[3])]

        ratio = Rmats._dpsi(ex1.ratio, ex1.ratio + ex2.ratio)

        return ratio

    @staticmethod
    def __A5SS_ratios(validation_dict, coords, gene_id, strand=None):
        ex1 = validation_dict['%s:%s-%s' % (gene_id, int(coords[0]) + 1, coords[1])]
        ex2 = validation_dict['%s:%s-%s' % (gene_id, int(coords[2]) + 1, coords[3])]

        ratio = Rmats._dpsi(ex1.ratio, ex1.ratio + ex2.ratio)

        return ratio

    @staticmethod
    def __RI_ratios(validation_dict, coords, gene_id, strand=None):
        pass

    @staticmethod
    def count_chg_genes(mats_dir, dpsi_thresh=0.2, pval_thresh=0.05):
        genes = set()
        mats_nn = 0
        for ev_type in Rmats.col_vals.keys():
            outputfile = '%s/%s.MATS.JC.txt' % (mats_dir, ev_type)
            #outputfile = '%s/%s.MATS.JCEC.txt' % (mats_dir, ev_type)
            with open(outputfile, 'r') as fp:
                for line in fp.readlines():
                    sline = line.split()
                    if sline[0] == "ID":
                        continue
                    else:
                        pvalue = float(sline[-5])
                        delta_psi = float(sline[-1])
                        comp = pvalue <= pval_thresh and abs(delta_psi) >= dpsi_thresh
                        if comp:
                            genes.add(sline[1].replace('"', ''))
        return genes

    @staticmethod
    def _count_sig_events(mats_dir, dpsi_thresh=0.2, pval_thresh=0.05, **kwargs):

        mats_nn = 0
        for ev_type in Rmats.col_vals.keys():
            outputfile = '%s/%s.MATS.JC.txt' % (mats_dir, ev_type)
            #outputfile = '%s/%s.MATS.JCEC.txt' % (mats_dir, ev_type)
            with open(outputfile, 'r') as fp:
                for line in fp.readlines():
                    sline = line.split()
                    if sline[0] == "ID":
                        continue
                    else:
                        evID = "%s:%s" %(ev_type, ":".join(sline[1:Rmats.col_vals[ev_type]['last_coord']+1]))
                        pvalue = float(sline[-5])
                        fdr = float(sline[-4])
                        delta_psi = float(sline[-1])
                        comp = pvalue <= pval_thresh and abs(delta_psi) >= dpsi_thresh
                        mats_nn += comp
        return mats_nn


    @staticmethod
    def _rank_rmats(mats_dir, dpsi_thresh=0.2, pval_thresh=0.05, step1_list=None, **kwargs):

        rank = []
        mats_nn = 0
        list_events = []
        ncounts = [0, 0, 0]
        only_pval = kwargs['only_pval']
        pval_f = kwargs['pval']
        ir = kwargs['ir']
        for ev_type in Rmats.col_vals.keys():
            if not ir and ev_type == 'RI':
                continue
            outputfile = '%s/%s.MATS.JC.txt' % (mats_dir, ev_type)
            #outputfile = '%s/%s.MATS.JCEC.txt' % (mats_dir, ev_type)
            with open(outputfile, 'r') as fp:
                for line in fp.readlines():
                    sline = line.split()
                    if sline[0] == "ID":
                        continue
                    else:
                        evID = "%s:%s" %(ev_type, ":".join(sline[1:Rmats.col_vals[ev_type]['last_coord']+1]))
                        pvalue = float(sline[-5])
                        fdr = float(sline[-4])
                        delta_psi = float(sline[-1])
                        # pass_thres = int(abs(delta_psi) >= dpsi_thresh and fdr <= pval_thresh)
                        if step1_list is None:
                            list_events.append(evID)
                        elif evID not in step1_list:
                            continue

                        flt1 = pvalue <= pval_thresh
                        flt2 = abs(delta_psi) >= dpsi_thresh
                        flt_all = pvalue <= pval_thresh and abs(delta_psi) >= dpsi_thresh

                        ncounts[0] += flt1
                        ncounts[1] += flt2
                        ncounts[2] += flt_all

                        comp = flt1 if only_pval else flt_all
                        mats_nn += comp
                        if comp:
                            rank.append([evID, delta_psi, pvalue])

        Rmats._print_stats(ncounts[0], ncounts[1], ncounts[2], dpsi_thresh=dpsi_thresh, pval_thresh=pval_thresh, method='rMats')
        rank = Rmats._sort_rank(rank, pval=pval_f)
        return rank, list_events

    """Public functions"""
    @staticmethod
    def rr_rank(dir1, dir2, dpsi_thresh=0.2, pval_thresh=0.05, **kwargs):
        rank1, ev_list = Rmats._rank_rmats(dir1, dpsi_thresh=dpsi_thresh, pval_thresh=pval_thresh, **kwargs)
        if len(rank1) == 0:
            print('Number of significant events in first comparision is 0. Exiting.')
            exit(-1)
        rank2, nlsvs2 = Rmats._rank_rmats(dir2, dpsi_thresh=dpsi_thresh, pval_thresh=pval_thresh,
                                          step1_list=ev_list, **kwargs)

        ratios = Rmats._rr_vals(rank1, rank2, method='rMats')
        return ratios

    @staticmethod
    def get_color(subsample=False):

        if subsample:
            color = Rmats.colors[1]
            lstyle = Rmats.lstyle_list[Rmats.plot_count % len(Rmats.lstyle_list)]
        else:
            color = Rmats.colors[Rmats.plot_count % len(Rmats.colors)]
            lstyle = '-'

        Rmats.plot_count += 1

        return color, lstyle

    @staticmethod
    def _validate_chg_genes(outputdir, validation_dict, values, **kwargs):
        
        events = Rmats._validate(outputdir, validation_dict, values, **kwargs)
        gn_dict = {}
        for ev, ev_vals in events.items():
            gn_id = ev.split(':')[1].replace('"','')
            try:
                gn_dict[gn_id].append(ev_vals)
            except:
                gn_dict[gn_id] = [ev_vals]
        return gn_dict


    @staticmethod
    def _validate(outputdir, validation_dict, values, **kwargs):
        print('VALIDATING rMATS files...')

        funcs = {'SE': Rmats.__SE_ratios,
                 'MXE': Rmats.__MXE_ratios,
                 'A5SS': Rmats.__A5SS_ratios,
                 'A3SS': Rmats.__A3SS_ratios,
                 'RI': Rmats.__RI_ratios}

        ev_dict = {}

        for ev_type in funcs.keys():
            if ev_type == 'RI':
                continue
            outputfile = '%s/%s.MATS.JC.txt' % (outputdir, ev_type)
            with open(outputfile, 'r') as fp:

                for line in fp.readlines():
                    if line.startswith('ID'):
                        continue
                    l = line.strip().split('\t')
                    evID = "%s:%s" %(ev_type, ":".join(l[1:Rmats.col_vals[ev_type]['last_coord']+1]))
                    evID = evID.replace('"', "")
#                    print(evID)
                    gene = l[1][1:-1]
                    strand = l[4]
                    delta = float(l[Rmats.col_vals[ev_type]['dPSI']])
                    pval = float(l[Rmats.col_vals[ev_type]['pval']])
                    FDR = float(l[Rmats.col_vals[ev_type]['FDR']])
                    
                    try:
                        ratio = funcs[ev_type](validation_dict[0], l[5:], gene)
                        ev_dict[evID] = (ratio, abs(delta), pval)

                    except KeyError:
                        values[Rmats.idx_stats['NOT_FOUND']] += 1
                        print ('Exon not found: %s Gene %s, junctions: %s' %(ev_type, gene, l[5:13]))
                        continue

        return ev_dict

    @staticmethod
    def _dpsi_delta(outputdir, validation_dict, dpsi_threshold=-1, **kwargs ):
        print('VALIDATING rMATS files...')

        funcs = {'SE': Rmats.__SE_ratios,
                 'MXE': Rmats.__MXE_ratios,
                 'A5SS': Rmats.__A5SS_ratios,
                 'A3SS': Rmats.__A3SS_ratios,
                 'RI': Rmats.__RI_ratios}

        deltas = []
        for ev_type in funcs.keys():
            if ev_type == 'RI':
                continue
            outputfile = '%s/%s.MATS.JC.txt' % (outputdir, ev_type)
            with open(outputfile, 'r') as fp:

                for line in fp.readlines():
                    if line.startswith('ID'):
                        continue
                    l = line.strip().split('\t')
                    evID = "%s:%s" % (ev_type, ":".join(l[1:Rmats.col_vals[ev_type]['last_coord'] + 1]))
                    gene = l[1][1:-1]
                    strand = l[4]
                    delta = abs(float(l[Rmats.col_vals[ev_type]['dPSI']]))
                    if dpsi_threshold != -1 and delta < dpsi_threshold: 
                        continue

                    FDR = float(l[Rmats.col_vals[ev_type]['FDR']])

                    try:
                        ratio = funcs[ev_type](validation_dict[0], l[5:], gene)
                        deltas.append(abs(delta-ratio))

                    except KeyError:
                        values[Rmats.idx_stats['NOT_FOUND']] += 1
                        print('Exon not found: %s Gene %s, junctions: %s' % (ev_type, gene, l[5:13]))
                        continue

        return deltas

    @staticmethod
    def rtpcr_dpsi(filename, list_of_pcr, output='./rtpcr_dpsi', **kwargs):
        data1 = []
        data2 = []

        vals_dict = {}
        labels = []
        with open(filename) as fp:
            for line in fp.readlines():
                if line.startswith('ID'):
                    continue
                tab = line.strip().split()

                inc_junc1 = "%s-%s" % (tab[8], int(tab[5]) + 1)
                exc_junc = "%s-%s" % (tab[8], int(tab[9]) + 1)
                inc_junc2 = "%s-%s" % (tab[6], int(tab[9]) + 1)

                eid1 = "%s:%s" % (tab[3], inc_junc1)
                eid2 = "%s:%s" % (tab[3], inc_junc2)
                eid3 = "%s:%s" % (tab[3], exc_junc)

                if eid1 in list_of_pcr[0]:
                    try:
                        vals_dict[eid1].append(float(tab[-1]))
                    except KeyError:
                        vals_dict[eid1] = [float(tab[-1])]

                elif eid2 in list_of_pcr[0]:
                    try:
                        vals_dict[eid2].append(float(tab[-1]))
                    except KeyError:
                        vals_dict[eid2] = [float(tab[-1])]

                # elif eid3 in list_of_pcr[1]:
                #     try:
                #         vals_dict[eid3].append(float(tab[-1]))
                #     except KeyError:
                #         vals_dict[eid3] = [float(tab[-1])]

        for eid, elst in vals_dict.items():
            labels.append(eid)
            data1.append(-float(list_of_pcr[0][eid][2]))
            data2.append(np.mean(elst))

        Rmats.dump_values_corr(data1, data2, output, labels=labels)
        print("N =", len(data1))
        Rmats.plot_corr(np.array(data1), np.array(data2), xax='RT-PCR DPSI', yax='rMATS DPSI',
                        title='%s N=%s' % ('rMATS', len(data1)), name=output)
        return

    @staticmethod
    def rtpcr_psi(filename, list_of_pcr, inv=False, output='./rtpcr_psi', title=''):
        return NotImplemented
