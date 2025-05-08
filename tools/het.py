from tools import Tools
import numpy as np
# from voila.api import Voila
import brewer2mpl
import csv

__author__    = "Jorge D. Vaquero"
__copyright__ = "Licensed under GNU GPLv3"

class Het(Tools):
    __pval_cache__ = {}

    plot_count = 0
    colors = np.array(brewer2mpl.get_map('Dark2', 'qualitative', 8).mpl_colors)[[0,1,2,3,4,5,6,7]]
    lstyle_list = ['solid', 'dashed', 'dotted', 'dashdot']
    fieldnames = ['gene_name ', 'gene_id',
                  'lsv_id', 'lsv_type',
                  'strand', 'seqid',
                  'Cerebellum_mean_psi', 'Muscle_mean_psi',
                  'INFOSCORE', 'TNOM',
                  'TTEST', 'WILCOXON',
                  'tnom_score', 'num_junctions',
                  'num_exons', 'de_novo_junctions',
                  'junctions_coords', 'exons_coords',
                  'ir_coords', 'ucsc_lsv_link']
    class Factory:
        def create(self): return Het()

    """private functions"""
    @staticmethod
    def __all_ratios(valdict, juncs, gene_id):

        ratios = np.array([valdict['%s:%s' % (gene_id, cc)].ratio for cc in juncs])
        sm = ratios.sum(axis=0)
        if sm[0] == 0 or sm[1] == 0:
            return None
        dpsi = [Het._dpsi(ratios[jidx], sm) for jidx in range(len(juncs))]
        return dpsi

    @staticmethod
    def __get_ratio(valdict, juncs, gene_id):
        ratios = np.array([valdict['%s:%s' % (gene_id, cc)].ratio for cc in juncs])
        #print(gene_id)
        sm = ratios.sum(axis=0)
        if sm[0] == 0 or sm[1] == 0:
            return None
        dpsi = 0
        for jidx in range(len(juncs)):
            dpsi = max(dpsi, Het._dpsi(ratios[jidx], sm))
        return dpsi

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
        dpsi = Het._dpsi(ratios[-1], sm)

        return dpsi

    @staticmethod
    def _matrix_area(matrix, V=0.2, absolute=True, collapsed_mat=False):
        """Returns the probability of an event to be above a certain threshold.
        The absolute flag describes if the value is absolute"""
        collapse = matrix if collapsed_mat else Het._collapse_matrix(matrix)
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
        collapsed = matrix if collapsed_mat else Het._collapse_matrix(matrix)
        xbins = np.linspace(-1, 1, num=collapsed.size+1)[:-1] + 1./collapsed.size
        if absolute: xbins = abs(xbins)
        return collapsed.dot(xbins)

    @staticmethod
    def _count_chg_genes(file1, dpsi_thresh=0.2, pval_thresh=0.05, statname=None):
        gene_set = set()
        with open(file1, 'r') as fp:
            reader = csv.DictReader(filter(lambda row: row[0]!='#', fp), delimiter='\t')
            psi_cols = [name for name in reader.fieldnames if 'median_psi' in name]
            if len(psi_cols) == 0:
                psi_cols = [name for name in reader.fieldnames if 'mean_psi' in name]

            for row in reader:
                lsv_id = row['lsv_id']
                gne_id = row['gene_id']
                gne_id = gne_id.replace('gene:', '')
                psi1 = [float(xx) for xx in row[psi_cols[0]].split(';')]
                psi2 = [float(xx) for xx in row[psi_cols[1]].split(';')]
                lsv_means = [abs(p2 - p1) >= dpsi_thresh for p1, p2 in zip(psi1, psi2)]
                area = [float(xx) <= pval_thresh for xx in row[statname].split(';')]
                valid = [x and y for (x, y) in zip(lsv_means, area)]
                if any(valid):
                    gene_set.add(gne_id)
        
        return gene_set


    @staticmethod
    def _count_sig_events(file1, dpsi_thresh=0.2, pval_thresh=0.05, statname=None, **kwargs):
        assert statname is not None
        lsv_rank = []
        majiq_nn = 0
        use_score = kwargs['use_score']
        use_overlap = kwargs['use_overlap']
        with open(file1, 'r') as fp:
            reader = csv.DictReader(filter(lambda row: row[0]!='#', fp), delimiter='\t')
            psi_cols = [name for name in reader.fieldnames if 'median_psi' in name]
            if len(psi_cols) == 0:
                psi_cols = [name for name in reader.fieldnames if 'mean_psi' in name]
            
    #    print(psi_cols)
            for row in reader:
                lsv_id = row['lsv_id']
                psi1 = [float(xx) for xx in row[psi_cols[0]].split(';')]
                psi2 = [float(xx) for xx in row[psi_cols[1]].split(';')]
                lsv_means = [abs(p2 - p1) for p1, p2 in zip(psi1, psi2)]  # vlsv.bins
                if len(lsv_means) > 2:
                    e_dpsi = [xx for xx in lsv_means]
                    junc_n = np.argmax(e_dpsi)
                else:
                    junc_n = 0

                if use_score:
                    colname = statname + '_score'
                    colname = colname.lower()
                else:
                    colname = statname
                area = [float(xx) for xx in row[colname].split(';')][junc_n]
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
                    lsv_exons = [[int(xx.split('-')[0]), int(xx.split('-')[1])] for xx in row['exons_coords'].split(';') if 'na' not in xx]
                    if np.any([ee in covered_exons for ee in lsv_exons if ee != lsv_exon_coords]):
                        continue
                    covered_exons.extend([ee for ee in lsv_exons if ee != lsv_exon_coords])
                if use_score:
                    bool_area = (area == 0)
                else:
                    bool_area = (area <= pval_thresh)
                mark_bool = abs(e_dpsi) >= dpsi_thresh and bool_area
                if mark_bool:
                    majiq_nn += 1

        return majiq_nn


    @staticmethod
    def _rank_majiq(vobj, dpsi_thresh=0.2, junc_selection={}, pval_thresh=0.05, step1_list=None, statname=None, **kwargs):

        assert statname is not None
        rank = []
        nlsv = -1
        total_events = []

#        for gid in vobj.get_gene_ids():
        lsv_rank = []
        nexp = 0
        nfdr = 0
        nexp_fdr = 0
        use_score = kwargs['use_score']
        pval = kwargs['pval']
        use_pval_med = kwargs['pval_median']
        use_overlap = kwargs['use_overlap']

        psi_cols = [name for name in vobj.fieldnames if 'median_psi' in name]
        if len(psi_cols) == 0:
            psi_cols = [name for name in reader.fieldnames if 'mean_psi' in name]
        for row in vobj:
            lsv_id = row['lsv_id']

            if step1_list is not None and lsv_id not in step1_list:
                continue
            nlsv += 1

            psi1 = [float(xx) for xx in row[psi_cols[0]].split(';')]
            psi2 = [float(xx) for xx in row[psi_cols[1]].split(';')]
            lsv_means = [p2 - p1 for p1, p2 in zip(psi1, psi2)]  # vlsv.bins
            if len(lsv_means) > 2:
                if step1_list is not None and lsv_id in junc_selection:
                    junc_n = junc_selection[lsv_id]
                else:
                    e_dpsi = [abs(xx) for xx in lsv_means]
                    junc_n = np.argmax(e_dpsi)
            else:
                junc_n = 0

            v_expected = lsv_means[junc_n]
            #area = 1.0 - Het._matrix_area(np.array(lsv_bins[junc_n]), dpsi_thresh, absolute=True, collapsed_mat=True)
            if use_score:
                colname = statname + '_score'
                colname = colname.lower()
            elif use_pval_med:
                colname = statname + '_quantile'
            else:
                colname = statname
            area = [float(xx) for xx in row[colname].split(';')][junc_n]
            #print(lsv_id, v_expected, area, junc_n)
            lsv_rank.append((row, v_expected, area, junc_n))
        lsv_rank.sort(key=lambda x: -x[2] if pval else abs(x[1]), reverse=True)

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

            if use_score:
                bool_area = (area == 0)
            else:
                bool_area = (area <= pval_thresh)

            mark_bool = abs(e_dpsi) >= dpsi_thresh and bool_area
            if mark_bool:
                rank.append(["%s#%d" % (lsv_id, junc_n), e_dpsi, area, bool_area])
            nexp += abs(e_dpsi) >= dpsi_thresh
            nfdr += bool_area
            nexp_fdr += mark_bool
            #print('[DEBUG] %s: %s (Delta median_psi = %.4f, PVAL = %.4f)' % (lsv_id, mark_bool, e_dpsi, area))

        Het._print_stats(nexp, nfdr, nexp_fdr, dpsi_thresh, pval_thresh, statname)
        rank = Het._sort_rank(rank, pval=pval)
        print ("Num of LSVs in majiq: %d" % nlsv)

        return rank, total_events

    """Public functions"""

    @staticmethod
    def rr_rank(file1, file2, dpsi_thresh=0.2, pval_thresh=0.05, statname=None, **kwargs):
        assert statname is not None

        with open(file1) as vobj:
        #with Voila(file1, 'r') as vobj:
            majiq_rank1, total_lsvs = Het._rank_majiq(csv.DictReader(filter(lambda row: row[0]!='#', vobj), delimiter='\t'), dpsi_thresh=dpsi_thresh,
                                                      pval_thresh=pval_thresh, statname=statname, **kwargs)
        if len(majiq_rank1) == 0:
            print('Number of significant events in first comparision is 0. Exiting.')
            exit(-1)

        junc_dict = {rr[0].split('#')[0]: int(rr[0].split('#')[1]) for rr in majiq_rank1}

        with open(file2) as vobj:
        # with Voila(file2, 'r') as vobj:
            majiq_rank2, total2 = Het._rank_majiq(csv.DictReader(filter(lambda row: row[0]!='#', vobj), delimiter='\t'), dpsi_thresh=dpsi_thresh,
                                                  pval_thresh=pval_thresh, junc_selection=junc_dict,
                                                  step1_list=total_lsvs, statname=statname, **kwargs)

        ratios = Het._rr_vals(majiq_rank1, majiq_rank2, method='Het')
        return ratios

    @staticmethod
    def get_color(subsample=False):

        if subsample:
            color = Het.colors[1]
            lstyle = Het.lstyle_list[Het.plot_count % len(Het.lstyle_list)]
        else:
            color = Het.colors[Het.plot_count % len(Het.colors)]
            lstyle = '-'

        Het.plot_count += 1
        return color, lstyle

    @staticmethod
    def _validate_chg_genes(outputfile, validation_dict, values, statname=None, **kwargs):
        events = Het._validate(outputfile, validation_dict, values, statname, **kwargs)
        gn_dict = {}
        for ev, ev_vals in events.items():
            gn_id = ':'.join(ev.split(':')[:-2])
            try:
                gn_dict[gn_id].append(ev_vals)
            except:
                gn_dict[gn_id] = [ev_vals]
        print("XX: ", gn_dict)
        return gn_dict

    @staticmethod
    def _validate(outputfile, validation_dict, values, statname=None, **kwargs):
        assert statname is not None
        print('VALIDATING MAJIQ files...')
        ev_dict = {}
        use_score = kwargs['use_score']
        ir = kwargs['ir']
        with open(outputfile, 'r') as fp:
            reader = csv.DictReader(filter(lambda row: row[0] != '#', fp), delimiter='\t')
            psi_cols = [name for name in reader.fieldnames if 'median_psi' in name]
            if len(psi_cols) == 0:
                psi_cols = [name for name in reader.fieldnames if 'mean_psi' in name]
            for row in reader:
                # if line.startswith('#'):
                #     continue

                # l = line.strip().split('\t')
                # print('KAKA', row)
                lsvtype = row['lsv_type']
                if ('i' in lsvtype[-1]) != ir:
                    continue
                # print('PE')
                gene = row['gene_id']
                juncs_l = row['junctions_coords'].split(';')
                psi1 = [float(xx) for xx in row[psi_cols[0]].split(';')]
                psi2 = [float(xx) for xx in row[psi_cols[1]].split(';')]
                dPSIs = [abs(p2 - p1) for p1, p2 in zip(psi1, psi2)]  # vlsv.bins
                if use_score:
                    colname = statname + '_score'
                    colname = colname.lower()
                else:
                    colname = statname
                probs = [float(x) for x in row[colname].split(';')]

                # psi1 = [abs(float(x)) for x in row['SRR1361714 median_psi'].split(';')]
                # psi2 = [abs(float(x)) for x in row['SRR1320579 median_psi'].split(';')]
                # dPSIs = [abs(x - psi2[xidx]) for xidx, x in enumerate(psi1)]


                njunc = len(dPSIs)
                if ir:
                    njunc -= 1

                juncs = [juncs_l[x] for x in range(njunc)]

                if ir:
                    lsv_chg = dPSIs[-1]
                    probs = probs[-1]
                    juncs.append(row['ir_coords'])

                else:
                    indx = np.argmax(dPSIs)
                    lsv_chg = dPSIs[indx]
                    pval = probs[indx]

                try:
                    if ir:
                        ratio = Het.__get_ratio_ir(validation_dict[1], juncs, gene)
                    else:
                        ratio = Het.__get_ratio(validation_dict[1], juncs, gene)

                    if ratio is None: continue
                    ev_dict[row['lsv_id']] = (ratio, abs(lsv_chg), pval)

                except KeyError as e:
                    values[Het.idx_stats['NOT_FOUND']] += 1
                    print ('Junction not found: Gene %s, junctions: %s %s' % (gene, juncs, e))

                    continue
        return ev_dict

    @staticmethod
    def _dpsi_delta(outputfile, validation_dict, dpsi_threshold=-1, statname=None, **kwargs):
        print('VALIDATING MAJIQ files...')
        ir = kwargs['ir']
        deltas = []
        with open(outputfile, 'r') as fp:
            reader = csv.DictReader(filter(lambda row: row[0] != '#', fp), delimiter='\t')
            psi_cols = [name for name in reader.fieldnames if 'median_psi' in name]
            if len(psi_cols) == 0:
                psi_cols = [name for name in reader.fieldnames if 'mean_psi' in name]
            for row in reader:
                lsvtype = row['lsv_type']
                if ('i' in lsvtype[-1]) != ir:
                    continue

                gene = row['gene_id']
                juncs_l = row['junctions_coords'].split(';')
                psi1 = [float(xx) for xx in row[psi_cols[0]].split(';')]
                psi2 = [float(xx) for xx in row[psi_cols[1]].split(';')]
                dPSIs = [abs(p2 - p1) for p1, p2 in zip(psi1, psi2)]

                try:
                    ratios = Het.__all_ratios(validation_dict[1], juncs_l, gene)
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
    def rtpcr_dpsi(filename, list_of_pcr, output='./rtpcr_dpsi', statname=None, **kwargs):
        data1 = []
        data2 = []

        dict_pcr = list_of_pcr[0]
        labels = []
        with open(filename, 'r') as fp:
            reader = csv.DictReader(filter(lambda row: row[0] != '#', fp), delimiter='\t')
            psi_cols = [name for name in reader.fieldnames if 'median_psi' in name]
            for row in reader:
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
                psi1 = [float(xx) for xx in row[psi_cols[0]].split(';')]
                psi2 = [float(xx) for xx in row[psi_cols[1]].split(';')]
                dPSIs = [(p2 - p1) for p1, p2 in zip(psi1, psi2)]  # vlsv.bins
                data2.append(float(dPSIs[junc_idx]))

        Het.dump_values_corr(data1, data2, output, labels=labels)
        print("N =", len(data1))
        Het.plot_corr(np.array(data1), np.array(data2), xax='RT-PCR DPSI', yax='%s DPSI' % statname,
                      title='%s N=%s' % (statname, len(data1)), name=output)

    @staticmethod
    def rtpcr_psi(filename, list_of_pcr, output='./rtpcr_psi', statname=None, **kwargs):

        data1_psi1 = []
        data2_psi1 = []
        data1_psi2 = []
        data2_psi2 = []

        dict_pcr = list_of_pcr[0]
        labels = []
        with open(filename, 'r') as fp:
            reader = csv.DictReader(filter(lambda row: row[0] != '#', fp), delimiter='\t')
            psi_cols = [name for name in reader.fieldnames if 'median_psi' in name]
            for row in reader:
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

                psi1 = [float(xx) for xx in row[psi_cols[0]].split(';')]
                data1_psi1.append(float(lsv[0]))
                data2_psi1.append(psi1[junc_idx])
        
                psi2 = [float(xx) for xx in row[psi_cols[1]].split(';')]
                data1_psi2.append(float(lsv[1]))
                data2_psi2.append(psi2[junc_idx])
        
                print("N =", len(data1_psi1))
                Het.plot_corr(np.array(data1_psi1), np.array(data2_psi1), xax='RT-PCR PSI1', yax='%s PSI1' % statname,
                              title='%s PSI2 N=%s' % (statname, len(data1_psi1)), name='./%s_psi2' % output)
        
                print("N =", len(data1_psi2))
                Het.plot_corr(np.array(data1_psi2), np.array(data2_psi2), xax='RT-PCR PSI2', yax='%s PSI2' % statname,
                              title='%s PSI2 N=%s' % (statname, len(data1_psi2)), name='./%s_psi2' % output)
        
    


