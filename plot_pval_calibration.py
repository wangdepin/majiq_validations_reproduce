import numpy as np
import argparse
import matplotlib.pyplot as plt
import csv

def cdf(nbins, divs):
    if len(divs) > 0:
        divs = np.array(divs)
        # Use the histogram function to bin the data
        counts, xbins = np.histogram(divs, bins=nbins, density=False)
        counts=counts.astype(float)/len(divs)
        # Now find the cdf
        cdf = np.cumsum(counts)
    return xbins, cdf


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('inputfile', nargs='+')
    parser.add_argument('-n', nargs='+')
    parser.add_argument('-o', required=True)
    parser.add_argument('--nbins', type=int, default=1000)
    args = parser.parse_args()

    xbins = np.linspace(0, 1, num=args.nbins + 1)
    
    stat_list = ['TNOM', 'WILCOXON', 'INFOSCORE', 'TTEST']

    nfiles = int(len(args.inputfile)/2)
    print( nfiles)
    fig, ax = plt.subplots(nfiles, 2, sharex=True, sharey=True)
    for i, fl in enumerate(args.inputfile):
        ev_list = {}
        with open(fl, 'r') as fp:
            reader = csv.DictReader(filter(lambda row: row[0]!='#', fp), delimiter='\t')
            psi_cols = [name for name in reader.fieldnames if 'mean_psi' in name]
            for row in reader:
                lsvtype = row['lsv_type']
#                if ('i' in lsvtype[-1]) != ir:
#                    continue
                # print('PE')
                gene = row['gene_id']
                if gene not in ev_list:
                    ev_list [gene] = {}
                    for statname in stat_list:
                        ev_list [gene][statname] = []

                    
                for statname in stat_list:
                    probs = [float(x) for x in row[statname].split(';')]

                    for p in probs:
                        ev_list[gene][statname].append(p)
   
        jn_idx = {}
        for gn, jnc_dict in ev_list.items():
            jn_idx[gn] = np.random.randint(0, len(jnc_dict['TNOM']))

        nr = int(i/2);
        nc = int(i%2);
        if nfiles == 1:
            f = ax[nc];
        else:
            f = ax[nr][nc];
        plots = []
        for idx, st in enumerate(stat_list):
            for gn, xx in ev_list.items():
                plots.append(ev_list[gn][st][jn_idx[gn]])
            

            f.plot(xbins[1:], cdf(args.nbins, plots)[1], label=st)
        f.plot(xbins,xbins, ls="--", c=".3")
        f.set_title(args.n[i])
        handles, labels = f.get_legend_handles_labels()
    plt.legend(handles, labels)

    plt.savefig(f'{args.o}', transparent=True)


if __name__ == '__main__':
    main()

