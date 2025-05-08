from tools import Tools
import brewer2mpl
import numpy as np
import gzip
from os import listdir
from os.path import isfile, join, basename

__author__    = "Jorge D. Vaquero"
__copyright__ = "Licensed under GNU GPLv3"

class Whippet(Tools):
    __pval_cache__ = {}

    colors = np.flip(np.array(brewer2mpl.get_map('RdPu', 'sequential', 9).mpl_colors))
    plot_count = 0
    lstyle_list = ['solid', 'dashed', 'dotted', 'dashdot']
    col_vals = {'pval': 8, 'dpsi': 7, 'coord': 2, 'node': 1, 'type': 4, 'gene': 0, 'strand': 3, 'CI_Width': 6,
                'inc_path': 11, 'exc_path': 12}
    col_vals_jnc = {'chrom': 0, 'start': 1, 'end': 2, 'nodes': 3, 'count': 4, 'strand': 5}
# Gene    Node    Coord   Strand  Type    Psi_A   Psi_B   DeltaPsi    Probability Complexity  Entropy
    class Factory:
        def create(self): return Whippet()

    """private functions"""
    @staticmethod
    def __get_id(row, jnc=False):
        return "%s:%s:%s" % (row[Whippet.col_vals['gene']], row[Whippet.col_vals['type']],
                             row[Whippet.col_vals['node']])

    @staticmethod
    def __get_psi_CI(whippet_dir):
        ci_pass = {}
        list_psi = ['%s/%s' %(whippet_dir, f) for f in listdir(whippet_dir)
                                                if isfile(join(whippet_dir, f)) and f.endswith('psi.gz')]
        for fn in list_psi:
            f = open(fn, 'rb')
            for line in gzip.open(f, 'rt').readlines()[1:]:
                sline = line.strip().split()
                ev_name = Whippet.__get_id(sline)
                try:
                    if float(sline[Whippet.col_vals['CI_Width']]) <= 0.1:
                        ci_pass[ev_name] += 1
                except KeyError:
                    ci_pass[ev_name] = 1
                except ValueError:
                    continue

        return ci_pass, len(list_psi)

    @staticmethod
    def __find_overlapping_exon(valdict, gene, start, end):
        for kk, vv in valdict[2][gene].items():
            s_ex, e_ex = (int(xx) for xx in kk.split('-'))
            if s_ex <= end and e_ex >= start:
                return vv
        raise KeyError

    @staticmethod
    def __extract_juncs(lst, gene_id, valdict):

        r = []
        for cc in lst[0]:
            k = '%s:%s' % (gene_id, cc)
            try:
                r.append(valdict[1][k].ratio)
            except KeyError as e:
                print('MISSING J: ', k)
                print(e)
                raise
        for cc in lst[1]:
            start, end = (int(xx) for xx in cc.split(':')[-1].split('-'))
            try:
                r.append(Whippet.__find_overlapping_exon(valdict, gene_id, start, end))
            except KeyError as e:
                print('MISSING EX: ', cc)
                raise
        return np.array(r)

    @staticmethod
    def count_chg_genes(whippet_dir, dpsi_thresh=0.2, pval_thresh=0.05, **kargs):

        genes = set()
        complx = True
        min_exp = 0.5

        whippet_dpsi = '%s/%s.diff.gz' % (whippet_dir, basename(whippet_dir))
        ci_pass, nfiles = Whippet.__get_psi_CI(whippet_dir)

        f = open(whippet_dpsi, 'rb')
        for line in gzip.open(f, 'rt').readlines()[1:]:
            sline = line.strip().split()

            dpsi = abs(float(sline[Whippet.col_vals['dpsi']]))
            pval = 1 - float(sline[Whippet.col_vals['pval']])
            ev_name = Whippet.__get_id(sline)
            pass_thres = pval <= pval_thresh and dpsi >= dpsi_thresh
            if ev_name not in ci_pass:
                continue

            complx_b = complx or sline[9] in ['K0', 'K1']
            if ci_pass[ev_name] >= (min_exp * nfiles) and complx_b and pass_thres:
               genes.add(sline[Whippet.col_vals['gene']])
        f.close()
        return genes


    @staticmethod
    def _count_sig_events(whippet_dir, dpsi_thresh=0.2, pval_thresh=0.05, **kargs):
        only_pval = kargs['onlypval']
        complx = True
        min_exp = kargs['min_exp']
        whippet_nn = 0

        whippet_dpsi = '%s/%s.diff.gz' % (whippet_dir, basename(whippet_dir))
        ci_pass, nfiles = Whippet.__get_psi_CI(whippet_dir)

        f = open(whippet_dpsi, 'rb')
        for line in gzip.open(f, 'rt').readlines()[1:]:
            sline = line.strip().split()

            dpsi = abs(float(sline[Whippet.col_vals['dpsi']]))
            pval = 1 - float(sline[Whippet.col_vals['pval']])
            ev_name = Whippet.__get_id(sline)
            pass_thres = pval <= pval_thresh
            if ev_name not in ci_pass:
                continue
            
            complx_b = complx or sline[9] in ['K0', 'K1']
            if ci_pass[ev_name] >= (min_exp *nfiles) and complx_b and pass_thres:
                whippet_nn += pass_thres
        f.close()
        return whippet_nn


    @staticmethod
    def _rank_whippet(whippet_dir, dpsi_thresh=0.2, pval_thresh=0.05, step1_list=None, **kwargs):

        rank = []
        whippet_nn = 0
        list_events = []
        ncounts = [0, 0, 0]
        pval_s = kwargs['pval']
        complx = kwargs['complex']
        min_exp = kwargs['min_exp']
        whippet_dpsi = '%s/%s.diff.gz' % (whippet_dir, basename(whippet_dir))

        ci_pass, nfiles = Whippet.__get_psi_CI(whippet_dir)

        f = open(whippet_dpsi, 'rb')
        for line in gzip.open(f, 'rt').readlines()[1:]:
            #print(line)
            sline = line.strip().split()

            dpsi = abs(float(sline[Whippet.col_vals['dpsi']]))
            pval = 1 - float(sline[Whippet.col_vals['pval']])
            ev_name = Whippet.__get_id(sline)
            pass_thres = int(abs(dpsi) >= dpsi_thresh)
            if ev_name not in ci_pass:
                continue
            if step1_list is None:
                list_events.append(ev_name)
            elif ev_name not in step1_list:
                continue
            
            complx_b = complx or sline[9] in ['K0', 'K1'];
            if ci_pass[ev_name] >= (min_exp *nfiles) and pass_thres and complx_b:
                whippet_nn += pass_thres
                rank.append([ev_name, dpsi, pval, pval <= pval_thresh])

            ncounts[0] += pval <= pval_thresh 
            ncounts[1] += dpsi >= dpsi_thresh
            ncounts[2] += pass_thres

        Whippet._print_stats(ncounts[0], ncounts[1], ncounts[2], dpsi_thresh=dpsi_thresh,
                             pval_thresh=pval_thresh, method='Whippet')
        rank = Whippet._sort_rank(rank, pval=pval)
        f.close()
        return rank, list_events

    """Public functions"""
    @staticmethod
    def rr_rank(file1, file2, dpsi_thresh=0.2, pval_thresh=0.05, **kwargs):

        rank1, event_list = Whippet._rank_whippet(file1, dpsi_thresh=dpsi_thresh, pval_thresh=pval_thresh, **kwargs)
        rank2, nlsvs2 = Whippet._rank_whippet(file2, dpsi_thresh=dpsi_thresh, pval_thresh=pval_thresh,
                                              step1_list=event_list, **kwargs)
        if len(rank1) == 0:
            print('Number of significant events in first comparision is 0. Exiting.')
            exit(-1)

        ratios = Whippet._rr_vals(rank1, rank2, method='Whippet')
        return ratios

    @staticmethod
    def get_color(subsample=False):

        if subsample:
            color = Whippet.colors[1]
            lstyle = Whippet.lstyle_list[Whippet.plot_count % len(Whippet.lstyle_list)]
        else:
            # color = Whippet.colors[(-1*Whippet.plot_count % len(Whippet.colors))]
            color = Whippet.colors[Whippet.plot_count % len(Whippet.lstyle_list)]
            lstyle = '-'
            print("COLOR: ", color)

        Whippet.plot_count += 1

        return color, lstyle

    @staticmethod
    def _validate_chg_genes(outputfile, validation_dict, values,  **kwargs):
        events = Whippet._validate(outputfile, validation_dict, values,  **kwargs)
        gn_dict = {}
        for ev, ev_vals in events.items():
            gn_id = ':'.join(ev.split(':')[:-2])
            try:
                gn_dict[gn_id].append(ev_vals)
            except:
                gn_dict[gn_id] = [ev_vals]
        return gn_dict

    @staticmethod
    def __parse_psi_path(sline, gene_id, j_dict, paths, valdict, exons, label):
        #print(gene_id, paths)
        for p in sline[Whippet.col_vals[label]].split(' '):
            #print('path: ', p)

            for p1 in p.split(','):
                temp_path = ([], set())
                xx = p1.split(':')[0].split('-')
                x1 = xx[0]
                for x2 in xx[1:]:
                    edge = "%s:%s-%s" % (gene_id, x1, x2)
                    #print(edge, end=' ')
                    if edge in j_dict:
                        c = "%s:%s" % (gene_id, j_dict[edge])
                        try:
                            temp_path[0].append(set(valdict[c].transcript_list))
                        except KeyError:
                            break
                    else:
                        temp_path[1].add("%s-%s" % (x1, x2))
                    x1 = x2
                else:

                    exons.append(temp_path[1])
                    paths.append(temp_path[0])
                    # if gene_id == 'ENSG00000160766':
                    #     print('EX: pre', exons)
                    #     print('PA: pre', paths)
                    #     print('TEMP: ', temp_path)
                    #     print('EX:', exons)
                    #     print('PA:', paths)


    @staticmethod
    def __parse_exonic_edges(path, exons, gene_d, nodes, valdict):

        for kk, vv in exons.items():
            gene_id = ':'.join(kk.split(':')[:-2])
            strand = gene_d[gene_id]
            for ii, pp in enumerate(vv):
                end = 0
                start = -1
                pre_c1 = [-1, -1]
                for xx in pp:
                    x1, x2 = xx.split('-')
                    ex1 = '%s:%s' % (gene_id, x1)
                    ex2 = '%s:%s' % (gene_id, x2)

                    try:
                        c1 = [int(cc) for cc in nodes[ex1].split('-')]
                    except KeyError:
                        if (start == -1):
                            continue
                        c1 = pre_c1

                    try:
                        c2 = [int(cc) for cc in nodes[ex2].split('-')]
                    except KeyError:
                        pre_c1 = c1
                        continue

                    if start == -1:
                        start = c1[0]
                        end = c1[1]

                    if (strand == '+' and c1[1] + 1 == c2[0]) or (strand == '-' and c2[1] + 1 == c1[0]):
                        start = min(min(c2[0], c1[0]), start)
                        end = max(max(c2[1], c1[1]), end)
                        continue
                    else:
                        if end != 0 and start != -1:
                            try:
                                v = Whippet.__find_overlapping_exon(valdict, gene_id, start, end)
                            except KeyError:
                                continue
                            path[kk][ii].append(set(v.transcript_list))
                        start = c2[0]
                        end = c2[1]

                if end != 0 and start != -1:
                    try:
                        v = Whippet.__find_overlapping_exon(valdict, gene_id, start, end)
                    except KeyError:
                        #path[kk][ii][0] = []
                        continue
                    path[kk][ii].append(set(v.transcript_list))

    @staticmethod
    def __get_ratio(valdict, j_inc, j_exc, gene_id):
        def __get_transcript_sum(trans_dict, paths, gene_id):
            sm = []
            if paths is not None:
                for pp in paths:
                    lp = pp.split('|')
                    for t in lp:
                        sm.append(trans_dict[-1][t])
            else:
                for t in trans_dict[3][gene_id]:
                    sm.append(trans_dict[-1][t])
            return np.array(sm)

        r_inc = __get_transcript_sum(valdict, j_inc, gene_id)
        if len(r_inc) == 0:
            raise KeyError
        sm_inc = r_inc.sum(axis=0)
        if sm_inc[0] == 0:
            raise KeyError

        r_exc = __get_transcript_sum(valdict, j_exc, gene_id)
        if len(r_exc) == 0:
           raise KeyError
        sm_exc = r_exc.sum(axis=0)
        dpsi = Whippet._dpsi(sm_inc, sm_inc + sm_exc)

        return dpsi

    @staticmethod
    def __merge_transcripts(in_p, out_p):
        for ev_id, plist in in_p.items():
            for pp in plist:
                if len(pp) == 0:
                    continue

                pp[0].intersection(*pp)
                try:
                    out_p[ev_id].add('|'.join(sorted(list(pp[0]))))
                except KeyError:
                    out_p[ev_id] = set()
                    out_p[ev_id].add('|'.join(sorted(list(pp[0]))))


    @staticmethod
    def _validate(whippet_dir, validation_dict, values, **kwargs):
        print('VALIDATING WHIPPET files...')
        ev_dict = {}
        complx = kwargs['complex']
        min_exp = kwargs['min_exp']
        whippet_dpsi = '%s/%s.diff.gz' % (whippet_dir, basename(whippet_dir))
        ir = kwargs['ir']

        list_junc = ['%s/%s' % (whippet_dir, f) for f in listdir(whippet_dir)
                     if isfile(join(whippet_dir, f)) and f.endswith('jnc.gz')]
        ci_pass = {}
        p_exc = {}
        p_inc = {}
        for fn in list_junc:
            print('PARSING %s', fn)
            j_dict = {}
            f = open(fn, 'rb')
            with gzip.open(f, 'rt') as gf:
                for line in gf.readlines()[1:]:
                    sline = line.strip().split()
                    junc_id = ':'.join(sline[Whippet.col_vals_jnc['nodes']].split(':')[:-1])
                    j_dict[junc_id] = '%s-%s' % (sline[Whippet.col_vals_jnc['start']], sline[Whippet.col_vals_jnc['end']])
            f.close()

            nodes = {}
            gene_d = {}
            ex_inc = {}
            ex_exc = {}
            t_exc = {}
            t_inc = {}
            fn_psi = fn[:-6] + 'psi.gz'
            f = open(fn_psi, 'rb')
            with gzip.open(f, 'rt') as gf:
                for line in gf.readlines()[1:]:
                    sline = line.strip().split()
                    ev_name = Whippet.__get_id(sline)
                    gene_id = sline[Whippet.col_vals['gene']]
                    gene_d[gene_id] = sline[Whippet.col_vals['strand']]

                    node_id = "%s:%s" % (sline[Whippet.col_vals['gene']], sline[Whippet.col_vals['node']])
                    nodes[node_id] = sline[Whippet.col_vals['coord']].split(':')[-1]

                    if ev_name not in ex_inc:
                        ex_inc[ev_name] = []
                        ex_exc[ev_name] = []
                        t_inc[ev_name] = []
                        t_exc[ev_name] = []

                    try:
                        if float(sline[Whippet.col_vals['CI_Width']]) <= 0.1:
                            ci_pass[ev_name] += 1
                        else:
                            continue
                    except KeyError:
                        ci_pass[ev_name] = 1
                    except ValueError:
                        continue

                    Whippet.__parse_psi_path(sline, gene_id, j_dict, t_inc[ev_name],
                                             validation_dict[1], ex_inc[ev_name], 'inc_path')
                    Whippet.__parse_psi_path(sline, gene_id, j_dict, t_exc[ev_name],
                                             validation_dict[1], ex_exc[ev_name], 'exc_path')

            Whippet.__parse_exonic_edges(t_inc, ex_inc, gene_d, nodes, validation_dict)
            Whippet.__parse_exonic_edges(t_exc, ex_exc, gene_d, nodes, validation_dict)
            f.close()

            Whippet.__merge_transcripts(t_inc, p_inc)
            Whippet.__merge_transcripts(t_exc, p_exc)

            del ex_inc
            del ex_exc
            del j_dict

        nfiles = len(list_junc)
        f = open(whippet_dpsi, 'rb')
        for line in gzip.open(f, 'rt').readlines()[1:]:
            sline = line.strip().split()

            dpsi = abs(float(sline[Whippet.col_vals['dpsi']]))
            pval = 1 - float(sline[Whippet.col_vals['pval']])
            ev_name = Whippet.__get_id(sline)
            if ev_name not in ci_pass or ci_pass[ev_name] < (min_exp * nfiles) or (sline[Whippet.col_vals['type']] =='RI' and not ir) or sline[Whippet.col_vals['type']] in ['BS', 'TS', 'TE']:
                continue

            complx_b = complx or sline[9] in ['K0', 'K1']
            gene_id = sline[Whippet.col_vals['gene']]
            node_id = '%s:%s' % (gene_id, sline[Whippet.col_vals['node']])

            try:
                if not complx_b or (ev_name not in p_exc and ev_name not in p_inc):
                    continue
                if ev_name in p_inc:
                    if ev_name not in p_exc:
                        inc = p_inc[ev_name]
                        exc = None
                    else:
                        inc = p_inc[ev_name]
                        exc = p_exc[ev_name]
                        #print(ev_name, p_exc[ev_name], p_inc[ev_name])
                        #print(line)
                    ratio = Whippet.__get_ratio(validation_dict, inc, exc, gene_id)
                ev_dict[ev_name] = (ratio, dpsi, pval)

            except KeyError as e:
                values[Whippet.idx_stats['NOT_FOUND']] += 1
                print('[%s] Junction not found: Gene %s' % (ev_name, gene_id))
                try:
                    print('Jinc: ', p_inc[ev_name])
                    print('Jexc: ', p_exc[ev_name])
                except KeyError:
                    continue
                print(e)
                continue

        f.close()
        return ev_dict

    @staticmethod
    def _dpsi_delta(whippet_dir, validation_dict, dpsi_threshold=-1, **kwargs):

        print('VALIDATING WHIPPET files...')
        deltas = []
        ir = kwargs['ir']
        complx = kwargs['complex']
        min_exp = kwargs['min_exp']
        complex = False
        whippet_dpsi = '%s/%s.diff.gz' % (whippet_dir, basename(whippet_dir))

        list_junc = ['%s/%s' % (whippet_dir, f) for f in listdir(whippet_dir)
                     if isfile(join(whippet_dir, f)) and f.endswith('jnc.gz')]
        ci_pass = {}
        p_exc = {}
        p_inc = {}
        for fn in list_junc:
            print('PARSING %s', fn)
            j_dict = {}
            f = open(fn, 'rb')
            with gzip.open(f, 'rt') as gf:
                for line in gf.readlines()[1:]:
                    sline = line.strip().split()
                    junc_id = ':'.join(sline[Whippet.col_vals_jnc['nodes']].split(':')[:-1])
                    j_dict[junc_id] = '%s-%s' % (sline[Whippet.col_vals_jnc['start']], sline[Whippet.col_vals_jnc['end']])
            f.close()

            nodes = {}
            gene_d = {}
            ex_inc = {}
            ex_exc = {}
            t_exc = {}
            t_inc = {}
            fn_psi = fn[:-6] + 'psi.gz'
            f = open(fn_psi, 'rb')
            with gzip.open(f, 'rt') as gf:
                for line in gf.readlines()[1:]:
                    sline = line.strip().split()
                    ev_name = Whippet.__get_id(sline)
                    gene_id = sline[Whippet.col_vals['gene']]
                    gene_d[gene_id] = sline[Whippet.col_vals['strand']]

                    node_id = "%s:%s" % (sline[Whippet.col_vals['gene']], sline[Whippet.col_vals['node']])
                    nodes[node_id] = sline[Whippet.col_vals['coord']].split(':')[-1]

                    if ev_name not in ex_inc:
                        ex_inc[ev_name] = []
                        ex_exc[ev_name] = []
                        t_inc[ev_name] = []
                        t_exc[ev_name] = []

                    try:
                        if float(sline[Whippet.col_vals['CI_Width']]) <= 0.1:
                            ci_pass[ev_name] += 1
                        else:
                            continue
                    except KeyError:
                        ci_pass[ev_name] = 1
                    except ValueError:
                        continue

                    Whippet.__parse_psi_path(sline, gene_id, j_dict, t_inc[ev_name],
                                             validation_dict[1], ex_inc[ev_name], 'inc_path')
                    Whippet.__parse_psi_path(sline, gene_id, j_dict, t_exc[ev_name],
                                             validation_dict[1], ex_exc[ev_name], 'exc_path')

            Whippet.__parse_exonic_edges(t_inc, ex_inc, gene_d, nodes, validation_dict)
            Whippet.__parse_exonic_edges(t_exc, ex_exc, gene_d, nodes, validation_dict)
            f.close()

            Whippet.__merge_transcripts(t_inc, p_inc)
            Whippet.__merge_transcripts(t_exc, p_exc)

            del ex_inc
            del ex_exc
            del j_dict

        nfiles = len(list_junc)
        f = open(whippet_dpsi, 'rb')

        for line in gzip.open(f, 'rt').readlines()[1:]:
            # print(line)
            sline = line.strip().split()

            dpsi = abs(float(sline[Whippet.col_vals['dpsi']]))
            if dpsi_threshold != -1 and dpsi < dpsi_threshold:
                continue
            ev_name = Whippet.__get_id(sline)
            if ev_name not in ci_pass or ci_pass[ev_name] < (min_exp * nfiles) or (sline[Whippet.col_vals['type']] =='RI' and not ir) or sline[Whippet.col_vals['type']] in ['BS', 'TS', 'TE']:
                continue

            complx_b = complx or sline[9] in ['K0', 'K1']
            gene_id = sline[Whippet.col_vals['gene']]
            node_id = '%s:%s' % (gene_id, sline[Whippet.col_vals['node']])

            try:
                if not complx_b or (ev_name not in p_exc and ev_name not in p_inc):
                    continue
                if ev_name in p_inc:
                    if ev_name not in p_exc:
                        inc = p_inc[ev_name]
                        exc = None
                    else:
                        inc = p_inc[ev_name]
                        exc = p_exc[ev_name]
                    ratio = Whippet.__get_ratio(validation_dict, inc, exc, gene_id)
                    if ratio == -1:
                        continue
                    deltas.append(abs(dpsi - ratio))

            except KeyError as e:
                print('[%s] Junction not found: Gene %s' % (ev_name, gene_id))
                try:
                    print('Jinc: ', p_inc[ev_name])
                    print('Jexc: ', p_exc[ev_name])
                except KeyError:
                    continue
                print(e)
                continue

        f.close()
        return deltas

    @staticmethod
    def _find_inc_node(inc_junc, exc_junc):
        print(inc_junc, exc_junc)
        t1 = inc_junc.split('-')
        t2 = exc_junc.split('-')

        if t1[0] == t2[0] and t1[1] != t2[1]:
            r = t1[1]
        elif t1[0] != t2[0] and t1[1] == t2[1]:
            r = t1[0]
        elif t1[0] != t2[0] and t1[1] != t2[1]:
            print('NOT SHARED NODE')
            r = 0
        return r



    @staticmethod
    def rtpcr_dpsi(whippet_dir, list_of_pcr, output='./rtpcr_dpsi', **kwargs):

        data1 = []
        data2 = []

        whippet_dpsi = '%s/%s.diff.gz' % (whippet_dir, basename(whippet_dir))

        list_junc = ['%s/%s' % (whippet_dir, f) for f in listdir(whippet_dir)
                     if isfile(join(whippet_dir, f)) and f.endswith('jnc.gz')]
        node_dict = {}
        j_dict = {}
        ev_dict = {}
        labels = []
        for fn in list_junc:
            found = {}
            f = open(fn, 'rb')
            with gzip.open(f, 'rt') as gf:
                for line in gf.readlines()[1:]:
                    sline = line.strip().split()
                    chrom = sline[Whippet.col_vals_jnc['chrom']]
                    junc_id = '%s:%s-%s' %(sline[Whippet.col_vals_jnc['chrom']], 
                                           sline[Whippet.col_vals_jnc['start']], 
                                           sline[Whippet.col_vals_jnc['end']])
                    gene_id = sline[Whippet.col_vals_jnc['nodes']].split(':')[0]
                    edge = sline[Whippet.col_vals_jnc['nodes']].split(':')[1]

                    if junc_id in list_of_pcr[0]:
                        exc_jid = "%s:%s" %(chrom, list_of_pcr[0][junc_id][3])
                        if exc_jid in found:
                            print(sline[Whippet.col_vals_jnc['nodes']])
                            n = Whippet._find_inc_node(edge, found[exc_jid])
                            node_id = "%s:%s" % (gene_id, n)
                            ev_dict[node_id] = (float(list_of_pcr[0][junc_id][2]), junc_id)
                        else:
                            found[junc_id] = edge
                    
                    elif junc_id in list_of_pcr[1]: 
                        inc_jid = "%s:%s" % (chrom, list_of_pcr[1][junc_id][3])
                        if inc_jid in found:
                            print(sline[Whippet.col_vals_jnc['nodes']])
                            n = Whippet._find_inc_node(found[inc_jid], edge)
                            node_id = "%s:%s" % (gene_id, n)
                            ev_dict[node_id] = (float(list_of_pcr[1][junc_id][2]), inc_jid)
                        else:
                            found[junc_id] = edge

        f = open(whippet_dpsi, 'rb')
        for line in gzip.open(f, 'rt').readlines()[1:]:
            # print(line)
            sline = line.strip().split()

            
            node_id = "%s:%s" %(sline[Whippet.col_vals['gene']], 
                                sline[Whippet.col_vals['node']])
            if node_id not in ev_dict:
                continue

            dpsi = -float(sline[Whippet.col_vals['dpsi']])
            labels.append(ev_dict[node_id][1])
            data1.append(ev_dict[node_id][0])
            data2.append(dpsi)

        print("N =", len(data1))
        Whippet.dump_values_corr(data1, data2, output, labels=labels)
        Whippet.plot_corr(np.array(data1), np.array(data2), xax='RT-PCR DPSI', yax='Whippet DPSI',
                        title='%s N=%s' % ('Whippet', len(data1)), name=output)
        return



