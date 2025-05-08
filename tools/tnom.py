from tools.het import Het


class Tnom(Het):
    class Factory:
        def create(self): return Tnom()

    @staticmethod
    def count_chg_genes(file1, dpsi_thresh=0.2, pval_thresh=0.05, statname=None):
        return Het._count_chg_genes(file1, dpsi_thresh=dpsi_thresh, pval_thresh=pval_thresh, statname='TNOM')

    @staticmethod
    def _count_sig_events(file1, dpsi_thresh=0.2, pval_thresh=0.05, **kwargs):
        return Het._count_sig_events(file1, dpsi_thresh=dpsi_thresh, pval_thresh=pval_thresh, statname='TNOM', **kwargs)

    @staticmethod
    def _rank_majiq(vobj, dpsi_thresh=0.2, junc_selection={}, pval_thresh=0.05, step1_list=None, **kwargs):
        return Het._rank_majiq(vobj, dpsi_thresh=dpsi_thresh, junc_selection=junc_selection, pval_thresh=pval_thresh,
                               step1_list=step1_list, statname='TNOM', **kwargs)
    @staticmethod
    def _validate_chg_genes(outputfile, validation_dict, vals, **kwargs):
        return Het._validate_chg_genes(outputfile, validation_dict, vals, statname='TNOM', **kwargs)

    @staticmethod
    def rr_rank(file1, file2, dpsi_thresh=0.2, pval_thresh=0.05, **kwargs):
        return Het.rr_rank(file1, file2, dpsi_thresh=dpsi_thresh, pval_thresh=pval_thresh, statname='TNOM', **kwargs)

    @staticmethod
    def _validate(outputfile, validation_dict, vals, **kwargs):
        return Het._validate(outputfile, validation_dict, vals, statname='TNOM', **kwargs)

    @staticmethod
    def _dpsi_delta(outputfile, validation_dict, dpsi_threshold=-1, **kwargs):
        return Het._dpsi_delta(outputfile, validation_dict, dpsi_threshold=-1, statname='TNOM', **kwargs)

    @staticmethod
    def rtpcr_dpsi(filename, list_of_pcr, output='./rtpcr_dpsi', **kwargs):
        return Het.rtpcr_dpsi(filename, list_of_pcr, output, statname='TNOM', **kwargs)

    @staticmethod
    def rtpcr_psi(filename, list_of_pcr, output='./rtpcr_psi', **kwargs):
        return Het.rtpcr_psi(filename, list_of_pcr, output, statname='TNOM', **kwargs)
