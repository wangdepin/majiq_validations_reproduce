import argparse
from tools import operator, all_stats
import matplotlib.pyplot as plt
import pandas
import upsetplot
import pickle
import os

__author__    = "Jorge D. Vaquero"
__copyright__ = "Licensed under GNU GPLv3"


def main(tools, pval=0.05, dpsi=0.2, title='Upset Plot', outname='upsetplot.pdf'):

    res = {}

    all_genes = set()
    tlist = list([x[1] for x in tools.keys()])
    for method, fil in tools.items():
        res[method[1]] = operator[method[0]].count_chg_genes(fil, dpsi_thresh=dpsi, pval_thresh=pval)


    rows = []

    sum_over = []
    #gene_list = {}
    for xx in range(len(tlist)):
        t1 = tlist[xx]
        # r = {x: True if t1 == x else False for x in tlist }
        # r['sum_over'] = len(res[t1])
        r = [t1]
        de = set()
        for mm in tlist:
            if mm != t1:
                de = de.union(res[mm])

        sg2 = res[t1].difference(de)
        if len(sg2) > 0:
            sum_over.append(len(sg2))
            rows.append(r)
            with open("%s_%s.tsv" %(outname, t1), 'w+') as of:
                [of.write("%s\n" % g) for g in sg2]

#            gene_list[t1] = sg2

        for i in range(xx+1, len(tlist)):
            lg = [t1]
            sg = res[t1]
            for yy in range(i, len(tlist)):
                #sg = sg.intersection(res[tlist[yy]])
                lg.append(tlist[yy])
                de = set()
                sg2 = res[t1]
                for mm in tlist:
                    if mm not in lg:
                        de = de.union(res[mm])
                    else:
                        sg2 = sg2.intersection(res[mm])
                sg2 = sg2.difference(de)
                if len(sg2) == 0:
                    continue
                rows.append(lg[:])
                # r = {x: True if x in lg else False for x in tlist}
                sum_over.append(len(sg2))
                # r['sum_over'] = len(sg)

    print(rows)
    print(sum_over)
    vals = upsetplot.from_memberships(rows, data=sum_over)

    print(vals)


    upsetplot.plot(vals, sort_by="cardinality")
    plt.title(title)
    current_figure = plt.gcf()

    current_figure.savefig("%s.pdf" % outname)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--tool', action='append', nargs=3, metavar=('tool_name', 'comp_id', 'file1'))
    parser.add_argument('-d', '--dpsi_threshold', action='store', dest='dpsi_t', default=0.2, type=float,
                        help='Threshold for minimum expected DPSI to be considered true changing')
    parser.add_argument('-p', '--pval_threshold', action='store', dest='pval_t', default=0.05, type=float,
                        help='Threshold for maximum pval to be considered true changing')
    parser.add_argument('--title', action='store', default='Upset Plot')
    parser.add_argument('-o', '--output', action='store', default='upsetplot')

    args = parser.parse_args()

    tools = {}
    print(all_stats)
    for xx in args.tool:
        if xx[0] not in all_stats:
            print('ERROR tool %s is not available' % xx[0])
            exit(-1)

        module_ = __import__('tools.' + xx[0].lower(), fromlist=xx[0].title())
        class_ = getattr(module_, xx[0].title())
        operator[xx[0]] = class_()
        tools[(xx[0], xx[1])] = xx[2]

    main(tools, pval=args.pval_t, dpsi=args.dpsi_t, title=args.title, outname=args.output)


