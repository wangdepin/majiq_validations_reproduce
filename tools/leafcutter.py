from tools import Tools
import numpy as np
import brewer2mpl
import os

__author__    = "Jorge D. Vaquero"
__copyright__ = "Licensed under GNU GPLv3"


class Leafcutter(Tools):
    __pval_cache__ = {}

    plot_count = 0

    """
    Colors and lstyle_list are example for global instance values to be used in get_color, depend on your 
    implementation of get_color, you don't really need them 
    """
    colors = np.array(brewer2mpl.get_map('Paired', 'qualitative', 11).mpl_colors)[[6, 7]]
    lstyle_list = ['solid', 'dashed', 'dotted', 'dashdot']

    class Factory:
        def create(self): return Leafcutter()

    """private functions"""
    @staticmethod
    def _rank_leafcutter(file_prefix, dpsi_thresh=0.2, pval_thresh=0.05, step1_list=None, **kwargs):

        clusters = {}
        rank = []
        lfctr_nn = 0
        list_events = []

        ncounts = [0, 0, 0]
        only_pval = kwargs['only_pval']
        pval_f = kwargs['pval']

        with open('%s_cluster_significance.txt' % file_prefix, 'r') as fp:
            for line in fp.readlines()[1:]:
                tab = line.strip().split()

                if tab[1] == 'Success' and float(tab[4]) <= pval_thresh:
                    clusters[tab[0]] = [float(tab[4]), 0, None]

        with open('%s_effect_sizes.txt' % file_prefix, 'r') as fp:
            for line in fp.readlines()[1:]:
                tab = line.strip().split()
                ids = tab[0].split(':')
                cluster_id = "%s:%s" % (ids[0], ids[3])
                dpsi = float(tab[4])
                if cluster_id not in clusters:
                    continue
                if abs(clusters[cluster_id][1]) < abs(dpsi):
                    clusters[cluster_id][1] = dpsi
                    clusters[cluster_id][2] = '%s-%s' % (ids[1], ids[2])

        for xx, vv in clusters.items():
            evID = "%s:%s" % (xx.split(':')[0], vv[2])
            #print (evID, step1_list)
            if step1_list is None:
                list_events.append(evID)
            elif evID not in step1_list:
                continue

            flt1 = vv[0] <= pval_thresh
            flt2 = abs(vv[1]) >= dpsi_thresh
            flt_all = vv[0] <= pval_thresh and abs(vv[1]) >= dpsi_thresh

            ncounts[0] += flt1
            ncounts[1] += flt2
            ncounts[2] += flt_all

            comp = flt1 if only_pval else flt_all
            lfctr_nn += comp
            if comp:
                rank.append([evID, vv[1], vv[0]])

        Leafcutter._print_stats(ncounts[0], ncounts[1], ncounts[2], dpsi_thresh=dpsi_thresh, pval_thresh=pval_thresh,
                                method='Leafcutter')
        rank = Leafcutter._sort_rank(rank, pval=pval_f)
        return rank, list_events

    @staticmethod
    def count_chg_genes(file_prefix, dpsi_thresh=0.2, pval_thresh=0.05):

        genes = set()
        genesymbol2ID = {}
        p = os.path.dirname(os.path.realpath(__file__))
        with open('%s/../data/geneid2name.tsv' % p) as fp:
            for ll in fp.readlines():
                tab = ll.strip().split()
                genesymbol2ID[tab[1]] = tab[0]

        clusters = {}
        with open('%s_cluster_significance.txt' % file_prefix, 'r') as fp:
            for line in fp.readlines()[1:]:
                tab = line.strip().split()
                if tab[1] == 'Success' and float(tab[4]) <= pval_thresh:
                    try:
                        gid = genesymbol2ID[tab[-1]]
                    except KeyError:
                        gid = tab[-1]
                    clusters[tab[0]] = [float(tab[4]), 0, gid]

        with open('%s_effect_sizes.txt' % file_prefix, 'r') as fp:
            for line in fp.readlines()[1:]:
                tab = line.strip().split()
                ids = tab[0].split(':')
                cluster_id = "%s:%s" % (ids[0], ids[3])
                dpsi = float(tab[4])
                if cluster_id not in clusters:
                    continue
                if abs(clusters[cluster_id][1]) < abs(dpsi):
                    clusters[cluster_id][1] = dpsi

        for xx, vv in clusters.items():
            if vv[0] <= pval_thresh and abs(vv[1]) >= dpsi_thresh:
                genes.add(vv[2])
        return genes

    """private functions"""
    @staticmethod
    def _count_sig_events(file_prefix, dpsi_thresh=0.2, pval_thresh=0.05, **kwargs):
        clusters = {}
        rank = []
        lfctr_nn = 0
        onlypval = kwargs['onlypval']
        with open('%s_cluster_significance.txt' % file_prefix, 'r') as fp:
            for line in fp.readlines()[1:]:
                tab = line.strip().split()
                if tab[1] == 'Success' and float(tab[4]) <= pval_thresh:
                    clusters[tab[0]] = [float(tab[4]), 0, None]

        with open('%s_effect_sizes.txt' % file_prefix, 'r') as fp:
            for line in fp.readlines()[1:]:
                tab = line.strip().split()
                ids = tab[0].split(':')
                cluster_id = "%s:%s" % (ids[0], ids[3])
                dpsi = float(tab[4])
                if cluster_id not in clusters:
                    continue
                if abs(clusters[cluster_id][1]) < abs(dpsi):
                    clusters[cluster_id][1] = dpsi
                    clusters[cluster_id][2] = '%s-%s' % (ids[1], ids[2])

        for xx, vv in clusters.items():
            if onlypval:
                comp = vv[0] <= pval_thresh
            else:
                comp = vv[0] <= pval_thresh and abs(vv[1]) >= dpsi_thresh

            lfctr_nn += comp

        return lfctr_nn

    @staticmethod
    def __all_ratios(valdict, juncs, chrom):
        ratios = np.array([valdict['%s:%s' % (chrom, cc[1])].ratio for cc in juncs])
        dpsi = [Leafcutter._dpsi(ratios[jidx], ratios.sum(axis=0)) for jidx in range(len(juncs))]
        return dpsi

    @staticmethod
    def __get_ratio(valdict, juncs, chrom):
        dpsi = 0
        ratios = np.array([valdict['%s:%s' % (chrom, cc[1])].ratio for cc in juncs])
        for jidx in range(len(juncs)):
            dpsi = max(dpsi, Leafcutter._dpsi(ratios[jidx], ratios.sum(axis=0)))

        return dpsi

    """Public functions"""
    @staticmethod
    def rr_rank(file_prefix1, file_prefix2, dpsi_thresh=0.2, pval_thresh=0.05, **kwargs):
        #pval=False, ir=False, only_pval=False, use_score=False):
        """
        rr_rank function Given 2 input files return the reproducibility ratio between the events in file1 and file2
        This reproducibility is based on the event name and type.
        :param file1: First deltapsi file used for RR
        :param file2: Second deltapsi file used for RR
        :param dpsi_thresh: E(DPSI) Threshold used for ranking
        :param pval_thresh: pval Threshold used for filtering (val<=pval_thresh )
        :param pval: Boolean to set if the rank is done based on pval sortint or E(DPSI) sorting. default false,
                     E[DPSI] sorting
        :return: ratio for reproducibility with N elements all passing filter of confidence or significance.
        """

        rank1, ev_list = Leafcutter._rank_leafcutter(file_prefix1, dpsi_thresh=dpsi_thresh, pval_thresh=pval_thresh,
                                                     **kwargs)
        if len(rank1) == 0:
            print('Number of significant events in first comparision is 0. Exiting.')
            exit(-1)

        rank2, nlsvs2 = Leafcutter._rank_leafcutter(file_prefix2, dpsi_thresh=dpsi_thresh, pval_thresh=pval_thresh,
                                                    step1_list=ev_list, **kwargs)

        ratios = Leafcutter._rr_vals(rank1, rank2, method='Leafcutter')
        return ratios

    @staticmethod
    def get_color(subsample=False):

        if subsample:
            color = Leafcutter.colors[1]
            lstyle = Leafcutter.lstyle_list[Leafcutter.plot_count % len(Leafcutter.lstyle_list)]
        else:
            color = Leafcutter.colors[Leafcutter.plot_count % len(Leafcutter.colors)]
            lstyle = '-'

            Leafcutter.plot_count += 1

        return color, lstyle

    @staticmethod
    def _validate_chg_genes(file_prefix, validation_dict, values, **kwargs):

        genesymbol2ID = {}
        clusters = {}
        p = os.path.dirname(os.path.realpath(__file__))
        with open('%s/../data/geneid2name.tsv' % p) as fp:
            for ll in fp.readlines():
                tab = ll.strip().split()
                genesymbol2ID[tab[1]] = tab[0]
        with open('%s_cluster_significance.txt' % file_prefix, 'r') as fp:
            for line in fp.readlines()[1:]:
                tab = line.strip().split('\t')
                if tab[4] == 'NA':
                    continue
                pval = float(tab[4])
                try:
                    clusters[tab[0]] = genesymbol2ID[tab[-1]]
                except:
                    glist = tab[-1].split(',')
                    for ii in range(len(glist)):
                        if glist[ii] in genesymbol2ID:
                            clusters[tab[0]] = genesymbol2ID[glist[ii]]
                            break;
                    else:
                        continue

        events = Leafcutter._validate(file_prefix, validation_dict, values, **kwargs)
        gn_dict = {}
        for ev, ev_vals in events.items():
            gn_id = clusters[ev]
            try:
                gn_dict[gn_id].append(ev_vals)
            except:
                gn_dict[gn_id] = [ev_vals]
        return gn_dict

    @staticmethod
    def _validate(file_prefix, validation_dict, values, **kwargs):
        print('VALIDATING Leafcutter files...')

        clusters = {}
        adapted_vali_dict = {}
        gene_to_chrom = {}
        with open('./gene_chrom.txt', 'r') as fp:
            for line in fp.readlines():
                tab = line.strip().split()
                gene_to_chrom[tab[1]] = tab[0]


        for xx,vv in validation_dict[1].items():
            tab = xx.split(':')
            id1 = ':'.join(tab[0:-1])
            new_id = '%s:%s' %(gene_to_chrom[id1], tab[-1])
            adapted_vali_dict[new_id] = vv

        with open('%s_cluster_significance.txt' % file_prefix, 'r') as fp:
            for line in fp.readlines()[1:]:
                tab = line.strip().split()
                if tab[1] == 'Success':
                    clusters[tab[0]] = [float(tab[4]), []]

        with open('%s_effect_sizes.txt' % file_prefix, 'r') as fp:
            for line in fp.readlines()[1:]:
                tab = line.strip().split()
                ids = tab[0].split(':')
                cluster_id = "%s:%s" % (ids[0], ids[3])
                dpsi = float(tab[4])
                if cluster_id not in clusters:
                    continue

                clusters[cluster_id][1].append([dpsi, '%s-%s' % (ids[1], ids[2])])
        ev_dict = {}

        for cluster_id, info in clusters.items():
            pval = info[0]
            juncs = info[1]

            delta = max([xx[0] for xx in juncs])
            
            try:
                ratio = Leafcutter.__get_ratio(adapted_vali_dict, juncs, cluster_id.split(':')[0])
                ev_dict[cluster_id] = (ratio, abs(delta), pval)

            except KeyError as e:
                print (e)
                values[Leafcutter.idx_stats['NOT_FOUND']] += 1
                print ('Junction not found: %s, junctions: %s' %(cluster_id, info))
                continue

        return ev_dict

    @staticmethod
    def _dpsi_delta(file_prefix, validation_dict, dpsi_threshold=-1, **kwargs):
        print('VALIDATING Leafcutter files...')

        clusters = {}
        adapted_vali_dict = {}
        gene_to_chrom = {}

        deltas = []

        with open('./gene_chrom.txt', 'r') as fp:
            for line in fp.readlines():
                tab = line.strip().split()
                gene_to_chrom[tab[1]] = tab[0]

        for xx, vv in validation_dict[1].items():
            tab = xx.split(':')
            id1 = ':'.join(tab[0:-1])
            new_id = '%s:%s' % (gene_to_chrom[id1], tab[-1])
            adapted_vali_dict[new_id] = vv

        with open('%s_cluster_significance.txt' % file_prefix, 'r') as fp:
            for line in fp.readlines()[1:]:
                tab = line.strip().split()
                if tab[1] == 'Success':
                    clusters[tab[0]] = [float(tab[4]), []]

        with open('%s_effect_sizes.txt' % file_prefix, 'r') as fp:
            for line in fp.readlines()[1:]:
                tab = line.strip().split()
                ids = tab[0].split(':')
                cluster_id = "%s:%s" % (ids[0], ids[3])
                dpsi = float(tab[4])
                if cluster_id not in clusters:
                    continue

                clusters[cluster_id][1].append([dpsi, '%s-%s' % (ids[1], ids[2])])
        ev_dict = {}

        for cluster_id, info in clusters.items():
            juncs = info[1]
            dlocal = -1

            try:
                ratios = Leafcutter.__all_ratios(adapted_vali_dict, juncs, cluster_id.split(':')[0])

            except KeyError as e:
                print(e)
                print('Junction not found: %s, junctions: %s' % (cluster_id, info))
                continue

            for indx, jnc in enumerate(juncs):
                delta = abs(jnc[0])
                if dpsi_threshold != -1 and delta < dpsi_threshold: 
                    continue
                d = abs(delta - ratios[indx])
                dlocal = max(d, dlocal)

            if dlocal >= 0:
                deltas.append(dlocal)

        return deltas


    @staticmethod
    def rtpcr_dpsi(filename, list_of_pcr, output='./rtpcr_dpsi', **kwargs):
        data1 = []
        data2 = []
        dict_pcr = list_of_pcr[0]
        #inv = kwargs['inv']
        inv = False
        vals_dict = {}
        labels = [] 
        with open(filename) as fp:
            for line in fp.readlines()[1:]:
                tab = line.strip().split()
                tt = tab[0].split(':')
                junc = "%s-%s" % (tt[1], tt[2])
                eid1 = "%s:%s" % (tt[0], junc)

                if eid1 in dict_pcr:
                    #print eid1, list_of_pcr[eid1][8], tab[-1]
                    try:
                        vals_dict[eid1].append(float(tab[-1]))
                    except KeyError:
                        vals_dict[eid1] = [float(tab[-1])]

        for eid, elst in vals_dict.items():
            labels.append(eid)
            data1.append(float(dict_pcr[eid][2]))
            data2.append(np.mean(elst))

        print("N =", len(data1))
        if inv:
            data2 = [-xx for xx in data2]
        Leafcutter.dump_values_corr(data1, data2, output, labels=labels)
        Leafcutter.plot_corr(np.array(data1), np.array(data2),
                             title='%s N=%s' % ('Leafcutter', len(data1)),
                             name=output)

    @staticmethod
    def rtpcr_psi(filename, list_of_pcr, output='./rtpcr_psi', **kwargs):

        dict_pcr = list_of_pcr[1]
        vals_dict_psi1 = {}
        vals_dict_psi2 = {}
        #inv = kwargs['inv']
        inv = False
        with open(filename) as fp:
            for line in fp.readlines()[1:]:
                tab = line.strip().split()
                tt = tab[0].split(':')
                junc = "%s-%s" % (tt[1], tt[2])
                eid1 = "%s:%s" % (tt[0], junc)

                if eid1 in dict_pcr:
                    try:
                        vals_dict_psi1[eid1].append(float(tab[2]))
                    except KeyError:
                        vals_dict_psi1[eid1] = [float(tab[2])]

                    try:
                        vals_dict_psi2[eid1].append(float(tab[3]))
                    except KeyError:
                        vals_dict_psi2[eid1] = [float(tab[3])]



        data1 = []
        data2 = []
        for eid, elst in vals_dict_psi1.items():
            data1.append(float(dict_pcr[eid][0]))
            data2.append(np.mean(elst))

        print ("PSI1 N =", len(data1))
        if inv:
            data2 = [-xx for xx in data2]
        Leafcutter.plot_corr(np.array(data1), np.array(data2), lims=(0, 1),
                             title='%s N=%s' % ('Leafcutter psi1', len(data1)),
                             name='./%s_psi1' % output)

        data1 = []
        data2 = []
        for eid, elst in vals_dict_psi2.items():
            data1.append(float(dict_pcr[eid][1]))
            data2.append(np.mean(elst))

        print ("PSI2 N =", len(data1))
        if inv:
            data2 = [-xx for xx in data2]
        Leafcutter.plot_corr(np.array(data1), np.array(data2), yax='Leafcutter', lims=(0, 1),
                             title='%s N=%s' % ('Leafcutter psi2', len(data1)),
                             name='./%s_psi2' % output)
