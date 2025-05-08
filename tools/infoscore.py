from tools.het import Het


class Infoscore(Het):
    class Factory:
        def create(self): return Infoscore()

    @staticmethod
    def count_chg_genes(file1, dpsi_thresh=0.2, pval_thresh=0.05, statname=None):
        return Het._count_chg_genes(file1, dpsi_thresh=dpsi_thresh, pval_thresh=pval_thresh, statname='INFOSCORE')

    @staticmethod
    def _count_sig_events(file1, dpsi_thresh=0.2, pval_thresh=0.05, **kwargs):
        return Het._count_sig_events(file1, dpsi_thresh=dpsi_thresh, pval_thresh=pval_thresh, statname='INFOSCORE', **kwargs)

    @staticmethod
    def _rank_majiq(vobj, dpsi_thresh=0.2, junc_selection={}, pval_thresh=0.05, step1_list=None, **kwargs):
        return Het._rank_majiq(vobj, dpsi_thresh=dpsi_thresh, junc_selection=junc_selection, pval_thresh=pval_thresh,
                               step1_list=step1_list, statname='INFOSCORE', **kwargs)
    @staticmethod
    def _validate_chg_genes(outputfile, validation_dict, values, **kwargs):
        return Het._validate_chg_genes(outputfile, validation_dict, values, statname='INFOSCORE', **kwargs)
    @staticmethod
    def rr_rank(file1, file2, dpsi_thresh=0.2, pval_thresh=0.05, **kwargs):
        return Het.rr_rank(file1, file2, dpsi_thresh=dpsi_thresh, pval_thresh=pval_thresh, statname='INFOSCORE', **kwargs)

    @staticmethod
    def _validate(outputfile, validation_dict, values, **kwargs):
        return Het._validate(outputfile, validation_dict, values, statname='INFOSCORE', **kwargs)

    @staticmethod
    def _dpsi_delta(files, validation_dict, dpsi_threshold=-1, **kwargs):
        return Het._dpsi_delta(files, validation_dict, dpsi_threshold=-1, statname='INFOSCORE', **kwargs)

    @staticmethod
    def rtpcr_dpsi(filename, list_of_pcr, output='./rtpcr_dpsi', **kwargs):
        return Het.rtpcr_dpsi(filename, list_of_pcr, output, statname='INFOSCORE', **kwargs)

    @staticmethod
    def rtpcr_psi(filename, list_of_pcr, output='./rtpcr_psi', **kwargs):
        return Het.rtpcr_psi(filename, list_of_pcr, output, statname='INFOSCORE', **kwargs)
