from tools import Tools
import numpy as np
import brewer2mpl

__author__    = "Jorge D. Vaquero"
__copyright__ = "Licensed under GNU GPLv3"

"""

This file is a Skeleton that can be used to create new tool files, please don't modify this one, just copy and 
modify the new one. 
The expected class name is the same as the file name if first character as upper case, you just need to fill up the 
public functions and search and replace the "Skeleton" name for your new class/file
There is no need of other changes through the code.

"""


class Skeleton(Tools):
    __pval_cache__ = {}

    plot_count = 0

    """
    Colors and lstyle_list are example for global instance values to be used in get_color, depend on your 
    implementation of get_color, you don't really need them 
    """
    colors = np.array(brewer2mpl.get_map('Paired', 'qualitative', 11).mpl_colors)[[0, 1]]
    lstyle_list = ['solid', 'dashed', 'dotted', 'dashdot']

    class Factory:
        def create(self): return Skeleton()

    """private functions"""

    #Any needed private functions

    """Public functions"""
    @staticmethod
    def rr_rank(file1, file2, dpsi_thresh=0.2, pval_thresh=0.05, pval=False):
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
        pass

    @staticmethod
    def get_color(subsample=False):
        """
        get_color returns the color to be used during different plots for this tool, usually using color global value
        :param subsample: In the Norton et al 2017 paper subsample is plot as different line type, this option change
        between color shades to fix color different line styles [optional]
        :return:
        """
        pass

    @staticmethod
    def validate(outputfile, validation_dict, dpsi_threshold=0.2, pval_threshold=0.05, permissive=False):
        """
        Validate matches the results from the simulated data runs for this tool.
        :param outputfile: tool output file or dict to be compared against truth
        :param validation_dict: 2 element list of validation object dictionary. Validation object can be found in
                                run_validation.py. First element of list is a dictionary of exonic Validation object
                                and its transcript inclusion. Second element are junction Validation objects, depends
                                of each tool, you can use one or another or both. Dictionary is hashed based on the
                                Validation id (gene, strand, coordinates)
        :param dpsi_threshold: E(DPSI) Threshold used to determine change or not change (val>=dpsi_threshold)
        :param pval_threshold: pval Threshold used for filtering (val<=pval_threshold )
        :param permissive: Allows tool gray values, where the truth is changing and the tool is not confident about
                           it but dpsi it is
        :return: List of up to 11 elements, each position meaning is a counter specified in Tools.idx_stats
         """
        pass
