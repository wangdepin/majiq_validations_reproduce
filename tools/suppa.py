from tools import Tools
import brewer2mpl
import numpy as np
import math
__author__    = "Jorge D. Vaquero"
__copyright__ = "Licensed under GNU GPLv3"

class Suppa(Tools):
    __pval_cache__ = {}

    plot_count = 0
    colors = np.array(brewer2mpl.get_map('Paired', 'qualitative', 11).mpl_colors)[[4, 5]]
    lstyle_list = ['solid', 'dashed', 'dotted', 'dashdot']
    col_vals = {'pval': 2, 'dpsi': 1, 'id': 0}

    class Factory:
        def create(self): return Suppa()

    """private functions"""
    @staticmethod
    def __SE_ratios(validation_dict, coords, gene_id, strand=None):

        e1 = coords[0].split('-')[0]
        s3 = coords[1].split('-')[1]

        c1a = validation_dict['%s:%s' % (gene_id, coords[0])]
        ac2 = validation_dict['%s:%s' % (gene_id, coords[1])]
        c1c2 = validation_dict['%s:%s-%s' % (gene_id, e1, s3)]

        ratio1 = Suppa._dpsi(c1a.ratio, c1a.ratio + c1c2.ratio)
        ratio2 = Suppa._dpsi(ac2.ratio, ac2.ratio + c1c2.ratio)

       #print 'SE', c1a.ratio, ac2.ratio, c1c2.ratio, ratio1, ratio2
        return max(ratio1, ratio2)

    @staticmethod
    def __MX_ratios(validation_dict, coords, gene_id, strand=None):
        c1a1 = validation_dict['%s:%s' % (gene_id, coords[0])]
        a1c2 = validation_dict['%s:%s' % (gene_id, coords[1])]
        c1a2 = validation_dict['%s:%s' % (gene_id, coords[2])]
        a2c2 = validation_dict['%s:%s' % (gene_id, coords[3])]

        ratio1 = Suppa._dpsi(c1a1.ratio, c1a1.ratio + c1a2.ratio)
        ratio2 = Suppa._dpsi(a2c2.ratio, a1c2.ratio + a2c2.ratio)

        return max(ratio1, ratio2)

    @staticmethod
    def __A5_ratios(validation_dict, coords, gene_id, strand=None):

        j1 = validation_dict['%s:%s' % (gene_id, coords[0])]
        j2 = validation_dict['%s:%s' % (gene_id, coords[1])]

        ratio = Suppa._dpsi(j1.ratio, j1.ratio + j2.ratio)

        return ratio

    @staticmethod
    def __A3_ratios(validation_dict, coords, gene_id, strand=None):
        j1 = validation_dict['%s:%s' % (gene_id, coords[0])]
        j2 = validation_dict['%s:%s' % (gene_id, coords[1])]

        ratio = Suppa._dpsi(j1.ratio, j1.ratio + j2.ratio)
        return ratio

    @staticmethod
    def __RI_ratios(validation_dict, coords, gene_id, strand=None):
        pass

    @staticmethod
    def __AF_ratios(validation_dict, coords, gene_id, strand=None):

        if strand == '+':
            idx1 = 1
            idx2 = 3
        else:
            idx1 = 0
            idx2 = 2

        j1 = validation_dict['%s:%s' % (gene_id, coords[idx1])]
        j2 = validation_dict['%s:%s' % (gene_id, coords[idx2])]

        ratio = Suppa._dpsi(j1.ratio, j1.ratio + j2.ratio)

        if strand == '-':
            ratio = 1 - ratio

        return ratio

    @staticmethod
    def __AL_ratios(validation_dict, coords, gene_id, strand=None):

        if strand == '+':
            idx1 = 0
            idx2 = 2
        else:
            idx1 = 1
            idx2 = 3

        j1 = validation_dict['%s:%s' % (gene_id, coords[idx1])]
        j2 = validation_dict['%s:%s' % (gene_id, coords[idx2])]

        ratio = Suppa._dpsi(j1.ratio, j1.ratio + j2.ratio)

        if strand == '+':
            ratio = 1 - ratio

        return ratio

    @staticmethod
    def count_chg_genes(suppa_file, dpsi_thresh=0.2, pval_thresh=0.05):

        genes = set()
        for line in open(suppa_file).readlines()[1:]:
            sline = line.strip().split()
            if sline[1] == 'nan':
                continue
            dpsi = float(sline[1])
            pval = float(sline[2])
            ev_name = sline[0]
            if pval <= pval_thresh and abs(dpsi) >= dpsi_thresh:
                gid = ev_name.split(';')[0]
                genes.add(gid)

        return genes


    @staticmethod
    def _count_sig_events(suppa_file, dpsi_thresh=0.2, pval_thresh=0.05, **kwargs):
        suppa_nn = 0
        onlypval = kwargs['onlypval']
        for line in open(suppa_file).readlines()[1:]:
            sline = line.strip().split()
            if sline[1] == 'nan':
                continue
            dpsi = float(sline[1])
            pval = float(sline[2])
            ev_name = sline[0]
            pass_thres = int(abs(dpsi) >= dpsi_thresh and pval <= pval_thresh)
            if onlypval:
                suppa_nn += (pval <= pval_thresh)
            else:
                suppa_nn += pass_thres
        return suppa_nn


    @staticmethod
    def _rank_suppa(suppa_file, dpsi_thresh=0.2, pval_thresh=0.05, step1_list=None, **kwargs):
        rank = []
        suppa_nn = 0
        list_events = []
        ncounts = [0, 0, 0]
        pval_f = kwargs['pval']
        for line in open(suppa_file).readlines()[1:]:
            sline = line.strip().split()
            if sline[1] == 'nan':
                continue
            dpsi = float(sline[1])
            pval = float(sline[2])
            ev_name = sline[0]
            pass_thres = int(abs(dpsi) >= dpsi_thresh and pval <= pval_thresh)
            if step1_list is None:
                list_events.append(ev_name)
            elif ev_name not in step1_list:
                continue
            if pval <= pval_thresh:
                suppa_nn += pass_thres
                rank.append([ev_name, dpsi, pval, pval <= pval_thresh])

            ncounts[0] += pval <= pval_thresh 
            ncounts[1] += abs(dpsi) >= dpsi_thresh
            ncounts[2] += pass_thres

        Suppa._print_stats(ncounts[0], ncounts[1], ncounts[2], dpsi_thresh=dpsi_thresh, pval_thresh=pval_thresh, method='Suppa')
        rank = Suppa._sort_rank(rank, pval=pval_f)
        return rank, list_events

    """Public functions"""
    @staticmethod
    def rr_rank(file1, file2, dpsi_thresh=0.2, pval_thresh=0.05, **kwargs):

        rank1, event_list = Suppa._rank_suppa(file1, dpsi_thresh=dpsi_thresh, pval_thresh=pval_thresh, **kwargs)
        rank2, nlsvs2 = Suppa._rank_suppa(file2, dpsi_thresh=dpsi_thresh, pval_thresh=pval_thresh,
                                          step1_list=event_list, **kwargs)
        if len(rank1) == 0:
            print('Number of significant events in first comparision is 0. Exiting.')
            exit(-1)

        ratios = Suppa._rr_vals(rank1, rank2, method='Suppa')
        return ratios

    @staticmethod
    def get_color(subsample=False):

        if subsample:
            color = Suppa.colors[1]
            lstyle = Suppa.lstyle_list[Suppa.plot_count % len(Suppa.lstyle_list)]
        else:
            color = Suppa.colors[Suppa.plot_count % len(Suppa.colors)]
            lstyle = '-'

        Suppa.plot_count += 1

        return color, lstyle

    @staticmethod
    def _validate_chg_genes(outputfile, validation_dict, values, **kwargs):

        events = Suppa._validate(outputfile, validation_dict, values, **kwargs)
        gn_dict = {}
        for ev, ev_vals in events.items():
            gn_id = ev.split(';')[0]
            try:
                gn_dict[gn_id].append(ev_vals)
            except:
                gn_dict[gn_id] = [ev_vals]
        return gn_dict

    @staticmethod
    def _validate(outputfile, validation_dict, values, **kwargs):
        print('VALIDATING SUPPA files...')

        funcs = {'SE': Suppa.__SE_ratios,
                 'MX': Suppa.__MX_ratios,
                 'A5': Suppa.__A5_ratios,
                 'A3': Suppa.__A3_ratios,
                 'RI': Suppa.__RI_ratios,
                 'AF': Suppa.__AF_ratios,
                 'AL': Suppa.__AL_ratios}


        #[Total Pos, Total Neg, TP, FP, TN, FN]
        ev_dict = {}
        with open(outputfile, 'r') as fp:

            for line in fp.readlines()[1:]:
                if line.startswith('Event_id'):
                    continue
                l = line.strip().split('\t')
                if l[1] == 'nan':  # assume this means not detected
                    continue
                pval = float(l[Suppa.col_vals['pval']])
                delta = float(l[Suppa.col_vals['dpsi']])
                gene_id, info = l[Suppa.col_vals['id']].split(';')

                mmmm = info.split(':')
                evtype = mmmm[0]
                chrom = mmmm[1]
                coords = mmmm[2:]

                if evtype == 'RI':
                    continue
                try:
                    ratio = funcs[evtype](validation_dict[1], coords[:-1], gene_id, strand=coords[-1])
                    ev_dict[l[Suppa.col_vals['id']]] = (ratio, abs(delta), pval)

                except KeyError as e:
                    values[Suppa.idx_stats['NOT_FOUND']] += 1
                    print ('[%s] Junction not found: Gene %s, junctions: %s' % (evtype, gene_id, coords[:-1]))
                    print (e)
                    continue
        return ev_dict

    @staticmethod
    def _dpsi_delta(outputfile, validation_dict, dpsi_threshold=-1, **kwargs):
        print('VALIDATING SUPPA files...')

        funcs = {'SE': Suppa.__SE_ratios,
                 'MX': Suppa.__MX_ratios,
                 'A5': Suppa.__A5_ratios,
                 'A3': Suppa.__A3_ratios,
                 'RI': Suppa.__RI_ratios,
                 'AF': Suppa.__AF_ratios,
                 'AL': Suppa.__AL_ratios}


        #[Total Pos, Total Neg, TP, FP, TN, FN]
        deltas = []
        with open(outputfile, 'r') as fp:

            for line in fp.readlines()[1:]:
                if line.startswith('Event_id'):
                    continue
                l = line.strip().split('\t')
                if l[1] == 'nan':  # assume this means not detected
                    continue
                delta = abs(float(l[Suppa.col_vals['dpsi']]))
                if dpsi_threshold != -1 and delta < dpsi_threshold: 
                    continue
                gene_id, info = l[Suppa.col_vals['id']].split(';')

                mmmm = info.split(':')
                evtype = mmmm[0]
                chrom = mmmm[1]
                coords = mmmm[2:]

                if evtype == 'RI':
                    continue
                try:
                    ratio = funcs[evtype](validation_dict[1], coords[:-1], gene_id, strand=coords[-1])
                    deltas.append(abs(delta - ratio))

                except KeyError as e:
                    print ('[%s] Junction not found: Gene %s, junctions: %s' % (evtype, gene_id, coords[:-1]))
                    print (e)
                    continue
        return deltas

    @staticmethod
    def rtpcr_dpsi(filename, list_of_pcr, output='./rtpcr_dpsi', **kwargs):
        data1 = []
        data2 = []
        dict_pcr = list_of_pcr[0]
        labels = []
        with open(filename) as fp:
            for line in fp.readlines()[1:]:
                if line.startswith('#'):
                    continue
                l = line.strip().split('\t')
                pval = float(l[Suppa.col_vals['pval']])
                delta = float(l[Suppa.col_vals['dpsi']])
                gene_id, info = l[Suppa.col_vals['id']].split(';')
                mmmm = info.split(':')
                evtype = mmmm[0]
                chrom = mmmm[1]
                coords = mmmm[2:]

                if evtype != 'SE' or math.isnan(delta): continue

                e1 = coords[0].split('-')[0]
                s3 = coords[1].split('-')[1]
                juncid1 = "%s:%s" % (chrom, coords[0])
                juncid2 = "%s:%s" % (chrom, coords[1])
                junc_exc = "%s-%s" % (e1, s3)

                if juncid1 in dict_pcr and dict_pcr[juncid1][3] == junc_exc:
                    labels.append(juncid1)
        #            print(juncid1, dict_pcr[juncid1], delta)
                    data1.append(float(dict_pcr[juncid1][2]))
                    data2.append(delta)
                elif juncid2 in dict_pcr and dict_pcr[juncid2][3] == junc_exc:
                    labels.append(juncid1)
         #           print(juncid2, dict_pcr[juncid2], delta)
                    data1.append(float(dict_pcr[juncid2][2]))
                    data2.append(delta)

        Suppa.dump_values_corr(data1, data2, output, labels=labels)
        print("N =", len(data1))
        Suppa.plot_corr(np.array(data1), np.array(data2), xax='RT-PCR DPSI', yax='SUPPA DPSI',
                        title='%s N=%s' % ('SUPPA', len(data1)), name=output)


    @staticmethod
    def rtpcr_psi(filename, list_of_pcr, inv=False, output='./rtpcr_psi', title=''):
        return NotImplemented
