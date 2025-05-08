from matplotlib import use
use('Agg')
import argparse
from tools import operator, all_stats
import csv
import sys
__author__    = "Jorge D. Vaquero"
__copyright__ = "Licensed under GNU GPLv3"

def parse_rtpcr_file(filename):
    list_of_pcr = [{}, {}]
    with open(filename, 'r') as fp:
        for row in csv.DictReader(filter(lambda row: row[0] != '#', fp), delimiter='\t'):
            key = "%s:%s" % (row['chromosome'], row['inc_junc'])
            list_of_pcr[0][key] = [row['rtpcr_psi1_avg'], row['rtpcr_psi2_avg'], row['rtpcr_dpsi'], row['exc_junc']]
            key = "%s:%s" % (row['chromosome'], row['exc_junc'])
            list_of_pcr[1][key] = [row['rtpcr_psi1_avg'], row['rtpcr_psi2_avg'], row['rtpcr_dpsi'], row['inc_junc']]
    return list_of_pcr

def main(rtpcr_file, tools, outFile, title='', ispsi=False):

    list_of_pcr = parse_rtpcr_file(rtpcr_file)
    for (method, mid), filename in tools.items():
        if ispsi:
            operator[method].rtpcr_psi(filename, list_of_pcr, output=outFile, title=title)
        else:
            operator[method].rtpcr_dpsi(filename, list_of_pcr, output=outFile, title=title)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('rtpcr_file', action='store', help='RTPCR values to compare with')
    parser.add_argument('-t', '--tool', action='append', nargs=3, metavar=('tool_name', 'id', 'file'))
    parser.add_argument('--title', action='store', dest='title')
    parser.add_argument('--psi',  action='store_true', default=False)

    parser.add_argument('-o', '--output', action='store', dest='outfile', default='./',
                        help='Output directory where the output files will be stored')

    args = parser.parse_args()
    print(" ".join(sys.argv))
    tools = {}
    for xx in args.tool:
        if xx[0] not in all_stats:
            print ('ERROR tool %s is not available' % xx[0])
            exit(-1)

        module_ = __import__('tools.' + xx[0].lower(), fromlist=xx[0].title())
        class_ = getattr(module_, xx[0].title())
        operator[xx[0]] = class_()
        tools[xx[0], xx[1]] = xx[2]

    main(args.rtpcr_file, tools, outFile=args.outfile, title=args.title, ispsi=args.psi)
