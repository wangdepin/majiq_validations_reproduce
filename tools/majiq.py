from tools import Tools
import numpy as np
# from voila.api import Voila
import brewer2mpl
import csv

__author__    = "Jorge D. Vaquero"
__copyright__ = "Licensed under GNU GPLv3"

class Majiq(Tools):
    __pval_cache__ = {}

    plot_count = 0
    colors = np.array(brewer2mpl.get_map('Paired', 'qualitative', 11).mpl_colors)[[0, 1, 8, 9]]
    lstyle_list = ['solid', 'dashed', 'dotted', 'dashdot']
    PSI_STRING = 'E(PSI) per LSV junction'

    class Factory:
        def create(self): return Majiq()

    """private functions"""
    @staticmethod
    def __all_ratios(valdict, juncs, gene_id):

        ratios = np.array([valdict['%s:%s' % (gene_id, cc)].ratio for cc in juncs])
        sm = ratios.sum(axis=0)
        if sm[0] == 0 or sm[1] == 0:
            return None
        dpsi = [Majiq._dpsi(ratios[jidx], sm) for jidx in range(len(juncs))]
        return dpsi

    @staticmethod
    def __get_ratio(valdict, juncs, gene_id):
        ratios = np.array([valdict['%s:%s' % (gene_id, cc)].ratio for cc in juncs])
        sm = ratios.sum(axis=0)
        if sm[0] == 0 or sm[1] == 0:
            return None
        dpsi = 0
        for jidx in range(len(juncs)):
            dpsi = max(dpsi, Majiq._dpsi(ratios[jidx], sm))
        return dpsi

    @staticmethod
    def __get_ratio_psis(valdict, juncs, gene_id):
        psis = []
        ratios = np.array([valdict['%s:%s' % (gene_id, cc)].ratio for cc in juncs])
        sm = ratios.sum(axis=0)
        if sm == 0:
            return -1
        for jidx in range(len(juncs)):
            psis.append(Majiq._psi(ratios[jidx], sm))
        return psis

    """private functions"""
    @staticmethod
    def __get_ratio_ir(valdict, juncs, gene_id):
        ircoord1, ircoord2 = [int(x) for x in juncs[-1].split('-')]
        irkey = 'IR:%s:%s-%s' % (gene_id, ircoord1+1, ircoord2-1)
        irkey = 'IR:%s:%s-%s' % (gene_id, ircoord1, ircoord2)
        if irkey in valdict:
            pass #        print(irkey, 'IN')
        else:
            print(irkey, 'intron not present')

        ratios = [valdict['%s:%s' % (gene_id, cc)].ratio for cc in juncs[:-1]]
        ratios.append(np.array(valdict[irkey].ratio))
        ratios = np.array(ratios)

        # with np.printoptions(precision=3, suppress=True):
        #     print(irkey, " ## 1 JJ FP", ratios[0], ratios[1])
        sm = ratios.sum(axis=0)
        if sm[0] == 0 or sm[1] == 0:
            return None
        dpsi = Majiq._psi(ratios[-1], sm)

        return dpsi

    @staticmethod
    def _matrix_area(matrix, V=0.2, absolute=True, collapsed_mat=False):
        """Returns the probability of an event to be above a certain threshold.
        The absolute flag describes if the value is absolute"""
        collapse = matrix if collapsed_mat else Majiq._collapse_matrix(matrix)
        # get the delta psi histogram borders based on the size of 'collapse'
        # grab the values inside the area of interest
        nbins = collapse.shape[0]
        delta_space = np.linspace(-1, 1, num=nbins + 1)[1:]
        if absolute:
            delta_space = abs(delta_space)
        if V < 0:
            border = abs(V) > delta_space
        else:
            border = V < delta_space
        area = collapse[border].sum()
        return area

    @staticmethod
    def _collapse_matrix(matrix):
        xbins, ybins = matrix.shape
        assert xbins == ybins
        DIAG = [matrix.diagonal(offset=xx).sum() for xx in range(1 - xbins, xbins)]
        return np.array(DIAG)

    @staticmethod
    def _expected_dpsi(matrix, collapsed_mat=False, absolute=True):
        """
        Calculate sum_dpsi=Prob(dpsi)*dpsi == sum_v = v*P(Delta PSI)
        """
        collapsed = matrix if collapsed_mat else Majiq._collapse_matrix(matrix)
        xbins = np.linspace(-1, 1, num=collapsed.size+1)[:-1] + 1./collapsed.size
        if absolute: xbins = abs(xbins)
        return collapsed.dot(xbins)

    @staticmethod
    def count_chg_genes(file1, dpsi_thresh=0.2, pval_thresh=0.05, **kwargs):
        gene_set = set()

        with open(file1, 'r') as fp:
            for row in csv.DictReader(filter(lambda row: row[0]!='#', fp), delimiter='\t'):
                lsv_id = row['lsv_id']
                gne_id = row['gene_id']
                gne_id = gne_id.replace('gene:', '')
                lsv_means = [abs(float(xx)) >= dpsi_thresh for xx in row['mean_dpsi_per_lsv_junction'].split(';')]
                area = [(1.0 - float(xx)) <= pval_thresh for xx in row['probability_changing'].split(';')]
                valid = [x and y for (x, y) in zip(lsv_means, area)]
                if any(valid):
                    gene_set.add(gne_id)
        return gene_set

    @staticmethod
    def _count_sig_events(file1, dpsi_thresh=0.2, pval_thresh=0.05, **kwargs):
        lsv_rank = []
        majiq_nn = 0
        use_overlap = kwargs['use_overlap']
        with open(file1, 'r') as fp:
            for row in csv.DictReader(filter(lambda row: row[0]!='#', fp), delimiter='\t'):
                lsv_id = row['lsv_id']
                lsv_means = [abs(float(xx)) for xx in row['mean_dpsi_per_lsv_junction'].split(';')]  # vlsv.bins
                if len(lsv_means) > 2:
                    e_dpsi = [abs(xx) for xx in lsv_means]
                    junc_n = np.argmax(e_dpsi)
                else:
                    junc_n = 0

                # area = 1.0 - Majiq._matrix_area(np.array(lsv_bins[junc_n]), dpsi_thresh, absolute=True, collapsed_mat=True)
                area = 1.0 - [float(xx) for xx in row['probability_changing'].split(';')][junc_n]
                if lsv_means[junc_n] <= dpsi_thresh: continue
                lsv_rank.append((row, lsv_means[junc_n], area, junc_n))

            lsv_rank.sort(key=lambda x: abs(x[1]), reverse=True)
            covered_exons = []
            for (row, e_dpsi, area, junc_n) in lsv_rank:
                lsv_id = row['lsv_id']
                if 'na' in lsv_id.split(':')[-1]:
                    continue
                if not use_overlap:
                    lsv_exon_coords = [int(coord) for coord in lsv_id.split(':')[-1].split('-')]
                    lsv_exons = [[int(xx.split('-')[0]), int(xx.split('-')[1])] for xx in row['exons_coords'].split(';') if 'na' not in xx ]
                    if np.any([ee in covered_exons for ee in lsv_exons if ee != lsv_exon_coords]):
                        continue
                    covered_exons.extend([ee for ee in lsv_exons if ee != lsv_exon_coords])
                mark_bool = abs(e_dpsi) >= dpsi_thresh and area <= pval_thresh
                if mark_bool:
                    majiq_nn += 1

        return majiq_nn


    @staticmethod
    def _rank_majiq(vobj, dpsi_thresh=0.2, junc_selection={}, pval_thresh=0.05, step1_list=None, **kwargs):

        rank = []
        nlsv = -1
        total_events = []

        lsv_rank = []
        nexp = 0
        nfdr = 0
        nexp_fdr = 0
        pval = kwargs['pval']
        use_overlap = kwargs['use_overlap']

        for row in vobj:
            lsv_id = row['lsv_id']

            if step1_list is not None and lsv_id not in step1_list:
                continue
            nlsv += 1

            lsv_means = [float(xx) for xx in row['mean_dpsi_per_lsv_junction'].split(';')] #vlsv.bins
            if len(lsv_means) > 2:
                if step1_list is not None and lsv_id in junc_selection:
                    junc_n = junc_selection[lsv_id]
                else:
                    e_dpsi = [abs(xx) for xx in lsv_means]
                    junc_n = np.argmax(e_dpsi)
            else:
                junc_n = 0

            v_expected = lsv_means[junc_n]
            #area = 1.0 - Majiq._matrix_area(np.array(lsv_bins[junc_n]), dpsi_thresh, absolute=True, collapsed_mat=True)
            area = 1.0 - [float(xx) for xx in row['probability_changing'].split(';')][junc_n]
            lsv_rank.append((row, v_expected, area, junc_n))

        lsv_rank.sort(key=lambda x: abs(x[1]), reverse=True)

        covered_exons = []
        for (row, e_dpsi, area, junc_n) in lsv_rank:
            lsv_id = row['lsv_id']
            if 'na' in lsv_id.split(':')[-1]:
                continue
            if not use_overlap:
                lsv_exon_coords = [int(coord) for coord in lsv_id.split(':')[-1].split('-')]
                lsv_exons = [[int(xx.split('-')[0]), int(xx.split('-')[1])] for xx in row['exons_coords'].split(';') if 'na' not in xx ]
                if np.any([ee in covered_exons for ee in lsv_exons if ee != lsv_exon_coords]):
                    continue
                covered_exons.extend([ee for ee in lsv_exons if ee != lsv_exon_coords])
            if step1_list is None:
                total_events.append(lsv_id)
            elif lsv_id not in step1_list:
                continue

            mark_bool = abs(e_dpsi) >= dpsi_thresh and area <= pval_thresh
            if mark_bool:
                rank.append(["%s#%d" % (lsv_id, junc_n), e_dpsi, area,  area <= pval_thresh])
            nexp += abs(e_dpsi) >= dpsi_thresh
            nfdr += area <= pval_thresh
            nexp_fdr += mark_bool

        Majiq._print_stats(nexp, nfdr, nexp_fdr, dpsi_thresh, pval_thresh, 'Majiq')
        rank = Majiq._sort_rank(rank, pval=pval)
        print ("Num of LSVs in majiq: %d" % nlsv)

        return rank, total_events

    """Public functions"""

    @staticmethod
    def rr_rank(file1, file2, dpsi_thresh=0.2, pval_thresh=0.05, **kargs):

        with open(file1) as vobj:
        #with Voila(file1, 'r') as vobj:
            majiq_rank1, total_lsvs = Majiq._rank_majiq(csv.DictReader(filter(lambda row: row[0]!='#', vobj), delimiter='\t'), dpsi_thresh=dpsi_thresh,
                                                        pval_thresh=pval_thresh, **kargs)
        if len(majiq_rank1) == 0:
            print('Number of significant events in first comparision is 0. Exiting.')
            exit(-1)

        junc_dict = {rr[0].split('#')[0]: int(rr[0].split('#')[1]) for rr in majiq_rank1}

        with open(file2) as vobj:
        # with Voila(file2, 'r') as vobj:
            majiq_rank2, total2 = Majiq._rank_majiq(csv.DictReader(filter(lambda row: row[0]!='#', vobj), delimiter='\t'), dpsi_thresh=dpsi_thresh,
                                                    pval_thresh=pval_thresh,
                                                    junc_selection=junc_dict,
                                                    step1_list=total_lsvs, **kargs)

        ratios = Majiq._rr_vals(majiq_rank1, majiq_rank2, method='Majiq')
        return ratios

    @staticmethod
    def get_color(subsample=False):

        if subsample:
            color = Majiq.colors[1]
            lstyle = Majiq.lstyle_list[Majiq.plot_count % len(Majiq.lstyle_list)]
        else:
            color = Majiq.colors[Majiq.plot_count % len(Majiq.colors)]
            lstyle = '-'

        Majiq.plot_count += 1
        return color, lstyle

    @staticmethod
    def _validate_chg_genes(outputfile, validation_dict, values, **kwargs):
        events = Majiq._validate(outputfile, validation_dict, values, **kwargs)
        gn_dict = {}
        for ev, ev_vals in events.items():
            gn_id = ':'.join(ev.split(':')[:-2])
            try:
                gn_dict[gn_id].append(ev_vals)
            except:
                gn_dict[gn_id] = [ev_vals]
        return gn_dict

    @staticmethod
    def _validate(outputfile, validation_dict, values, **kwargs):
        print('VALIDATING MAJIQ files...')
        ir = kwargs['ir']
        ev_dict = {}
        with open(outputfile, 'r') as fp:
            for row in csv.DictReader(filter(lambda row: row[0]!='#', fp), delimiter='\t'):

                lsvtype = row['lsv_type']
                if ('i' in lsvtype[-1]) != ir:
                    continue
                gene = row['gene_id']
                juncs_l = row['junctions_coords'].split(';')
                dPSIs = [abs(float(x)) for x in row['mean_dpsi_per_lsv_junction'].split(';')]
                probs = [float(x) for x in row['probability_changing'].split(';')]
                njunc = len(dPSIs)
                if ir:
                    njunc -= 1

                juncs = [juncs_l[x] for x in range(njunc)]

                if ir:
                    lsv_chg = dPSIs[-1]
                    probs = 1 - probs[-1]
                    juncs.append(row['ir_coords'])
                    #if(len(juncs)) > 2 : continue

                else:
                    indx = np.argmax(dPSIs)
                    lsv_chg = dPSIs[indx]
                    probs = 1 - probs[indx]

                try:
                    if ir:
                        ratio = Majiq.__get_ratio_ir(validation_dict[1], juncs, gene)
                    else:
                        ratio = Majiq.__get_ratio(validation_dict[1], juncs, gene)

                    if ratio is None: continue
                    
                    ev_dict[row['lsv_id']] = (ratio, abs(lsv_chg), probs)

                except KeyError as e:
                    values[Majiq.idx_stats['NOT_FOUND']] += 1
                    print ('Junction not found: Gene %s, junctions: %s %s' %(gene, juncs, e))

                    continue

        return ev_dict

    @staticmethod
    def _dpsi_delta(outputfile, validation_dict, dpsi_threshold=-1, **kwargs ):
        print('VALIDATING MAJIQ files...')
        ir = kwargs['ir']
        deltas = []
        with open(outputfile, 'r') as fp:
            for row in csv.DictReader(filter(lambda row: row[0] != '#', fp), delimiter='\t'):
                lsvtype = row['lsv_type']
                if ('i' in lsvtype[-1]) != ir:
                    continue
                gene = row['gene_id']
                juncs_l = row['junctions_coords'].split(';')
                dPSIs = [abs(float(x)) for x in row['mean_dpsi_per_lsv_junction'].split(';')]
                probs = [float(x) for x in row['probability_changing'].split(';')]

                try:
                    ratios = Majiq.__all_ratios(validation_dict[1], juncs_l, gene)
                    if ratios is None: continue
                except KeyError as e:
                    print('Junction not found: Gene %s, junctions: %s %s' % (gene, juncs_l, e))
                    continue

                delta_local = -1
                for indx, junc in enumerate(juncs_l):
                    lsv_chg = abs(dPSIs[indx])
                    if dpsi_threshold != -1 and lsv_chg < dpsi_threshold:
                        continue
                    delta_local = max(abs(lsv_chg - ratios[indx]), delta_local)

                if delta_local >= 0: 
                    deltas.append(delta_local)
        return deltas

    @staticmethod
    def translate_quantifications(infile, validation_dict, outfile):
        print('Translating MAJIQ files...')
        values = [0] * len(Majiq.idx_stats)
        with open(infile, 'r') as fp:
            reader = csv.DictReader(filter(lambda row: row[0]!='#', fp), delimiter='\t')
            with open(outfile, "w") as f_out:
                writer = csv.DictWriter(f_out, fieldnames=reader.fieldnames, delimiter="\t")
                writer.writeheader()
                for row in reader:
                    gene = row['gene_id']
                    juncs_l = row['junctions_coords'].split(';')
                    PSIs = [abs(float(x)) for x in row[Majiq.PSI_STRING].split(';')]

                    njunc = len(PSIs)
                    juncs = [juncs_l[x] for x in range(njunc)]

                    try:
                        ratio = Majiq.__get_ratio_psis(validation_dict[1], juncs, gene)
                        if ratio is None:
                            continue
                        row[Majiq.PSI_STRING] = ';'.join(['%.3f' % xx for xx in ratio])
                        writer.writerow(row)

                    except KeyError as e:
                        values[Majiq.idx_stats['NOT_FOUND']] += 1
                        print('Junction not found: Gene %s, junctions: %s %s' % (gene, juncs, e))

                        continue
        return

    @staticmethod
    def rtpcr_dpsi(filename, list_of_pcr, output='./rtpcr_dpsi', **kwargs):

        data1 = []
        data2 = []

        dict_pcr = list_of_pcr[0]
        labels = []
        with open(filename, 'r') as fp:
            for row in csv.DictReader(filter(lambda row: row[0] != '#', fp), delimiter='\t'):
                list_juncs = row['junctions_coords'].split(';')
                for dx, jj in enumerate(list_juncs):
                    eidx = '%s:%s' % (row['seqid'], jj)
                    if eidx in dict_pcr:
                        lsv = dict_pcr[eidx]
                        junc_idx = dx
                        if dict_pcr[eidx][3] in list_juncs:
                            break
                else:
                    continue

                labels.append(eidx)
                data1.append(float(lsv[2]))
                data2.append(float(row['mean_dpsi_per_lsv_junction'].split(';')[junc_idx]))
        Majiq.dump_values_corr(data1, data2, output, labels=labels)

        print("N =", len(data1))
        Majiq.plot_corr(np.array(data1), np.array(data2), xax='RT-PCR DPSI', yax='MAJIQ DPSI',
                        title='%s N=%s' % ('MAJIQ', len(data1)), name=output)

    @staticmethod
    def rtpcr_psi(filename, list_of_pcr, output='./rtpcr_psi', **kwargs):

        data1_psi1 = []
        data2_psi1 = []

        data1_psi2 = []
        data2_psi2 = []

        dict_pcr = list_of_pcr[0]
        labels = []
        with open(filename, 'r') as fp:
            for row in csv.DictReader(filter(lambda row: row[0] != '#', fp), delimiter='\t'):
                junc_idx = None
                list_juncs = row['junctions_coords'].split(';')
                for dx, jj in enumerate(list_juncs):
                    eidx = '%s:%s' % (row['seqid'], jj)
                    if eidx in dict_pcr:
                        lsv = dict_pcr[eidx]
                        junc_idx = dx
                        if dict_pcr[eidx][3] in list_juncs:
                            break
                else:
                    continue
                
                labels.append(eidx)
                psi_cols = [name for name in row.fieldnames if 'mean_psi' in name]
                psi1 = [float(xx) for xx in row[psi_cols[0]].split(';')]
                data1_psi1.append(float(lsv[0]))
                data2_psi1.append(psi1[junc_idx])

                psi2 = [float(xx) for xx in row[psi_cols[1]].split(';')]
                data1_psi2.append(float(lsv[1]))
                data2_psi2.append(psi2[junc_idx])
        
        Majiq.dump_values_corr(data1, data2, output, labels=labels)
        
        print("N =", len(data1_psi1))
        Majiq.plot_corr(np.array(data1_psi1), np.array(data2_psi1), xax='RT-PCR PSI1', yax='MAJIQ PSI1',
                        title='MAJIQ PSI1 N=%s' %  len(data1_psi1), name='./%s_psi2' % output)

        print("N =", len(data1_psi2))
        Majiq.plot_corr(np.array(data1_psi2), np.array(data2_psi2), xax='RT-PCR PSI2', yax='MAJIQ PSI2',
                        title='MAJIQ PSI2 N=%s' % len(data1_psi2), name='./%s_psi2' % output)




