import argparse
import matplotlib.pyplot as plt
import numpy as np
import csv
import random
import subprocess
import os
from collections import defaultdict
from enum import Enum

def load_junction_filter(infile):
    output = defaultdict(list)
    with open(infile) as fp:
        reader = csv.DictReader(fp, dialect='excel-tab')
        for row in reader:
            output[row['Gene ID']].append([int(x) for x in row['Junctions coords'].split('-')])
    return output

def pdf_or_cdf(xbins, divs):
    if len(divs) > 0:
        divs = np.array(divs)
        print(divs)
        distribs = np.less.outer(divs, xbins).mean(axis=(0, 1))
        print(distribs)
        return distribs
#        return len(divs), distribs

def cdf(nbins, divs):
    if len(divs) > 0:
#        divs = np.array(divs).mean(axis=1)
        divs = np.array(divs).max(axis=1)
        # Use the histogram function to bin the data
        counts, xbins = np.histogram(divs, bins=nbins, density=False)
        counts=counts.astype(float)/len(divs)
        # Now find the cdf
        cdf = np.cumsum(counts)
    return xbins, cdf

def epsi_complexity_filter(psi_tsv, k, cmplx, last_one):
    def filtfn(psi_mat):
        psi_mat = np.array(psi_mat)
        nexp, njunc = psi_mat.shape
        return nexp == k #and (cmplx is None or (njunc >= cmplx if last_one else njunc == cmplx))

    for row in filter(filtfn, psi_tsv.values()):
        e_psi = np.array(row)
        divs = abs(e_psi - np.median(e_psi, axis=0))
        best_i = divs.sum(axis=0).argmax()
        yield divs[:, best_i]

def voila_epsi_reader(voila_files, gene_ids=None, nfiles=4):
    psi = defaultdict(lambda: defaultdict(list))
    seen_genes = {}
    for voila_file in random.sample(voila_files, nfiles):
        with open(voila_file) as fp:
            print('voila_file', voila_file)
            reader = csv.DictReader(fp, dialect='excel-tab')
            for row in reader:
                gene_id = row['Gene ID']
                lsv_id = row['LSV ID']
                psi_vals = list(map(float, row['E(PSI) per LSV junction'].split(';')))
                psi[gene_id][lsv_id].append(psi_vals)

    gene_psi = [] 
    for gn, pslist in psi.items():
        mx = 0.0
        val = None
        for lsv, data in pslist.items():
            if len(data)< nfiles :
                continue
            e_psi = np.array(data)
   #         print('####')
   #         print(e_psi)
            divs = abs(e_psi - np.median(e_psi, axis=0))
   #         print(divs)
            best_i = divs.sum(axis=0).argmax()
   #         print(best_i)
            tmp = divs[:, best_i].max()
            
            if tmp > mx:
                val = divs[:, best_i]
                mx = tmp
        if val is None:
            continue
        gene_psi.append(val)

            
    print('n lsv:', len(gene_psi))
    return gene_psi


def get_psi_tsv(psi_tsv: np.ndarray, k, complexity=(None,)):
    for cmplx in complexity:
        yield list(epsi_complexity_filter(psi_tsv, k, cmplx, cmplx == complexity[-1]))

def get_epsi_dists(voila_files, args, xbins, gene_ids=None, name='out.npz'):
    psi = voila_epsi_reader(voila_files, gene_ids=gene_ids, nfiles=args.nrepl)
    np.save(name, np.array(psi))
    xbins, res = cdf(args.nbins, psi)
    return res


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('-t', '--tool', action='append', nargs='+', help='list of curves to be ploted, first string in -t is the identifier')
    parser.add_argument('-o', required=True)
    parser.add_argument('-n', type=int, default=1)
    parser.add_argument('-k', dest="nrepl", type=int, default=4)
    parser.add_argument('--nbins', type=int, default=1000)
    parser.add_argument('--nthreads', type=int, default=16)
    parser.add_argument('--complexity', type=int, nargs='+', default=(None,))
    parser.add_argument('--gene-ids-file', default=None)
    args = parser.parse_args()

    xbins = np.linspace(0, 1, num=args.nbins + 1)
    gene_ids=None
    if args.gene_ids_file is not None:
        gene_ids = list(line.strip() for line in open(args.gene_ids_file))

    for xx in args.tool:
        distribs =[ get_epsi_dists(xx[1:], args, xbins, gene_ids, name=xx[0]) for it in range(args.n)]

        distrib_array = np.array(distribs)
        #mu = distrib_array.mean(axis=0)
        mu = distrib_array.max(axis=0)
        std = distrib_array.std(axis=0)
        lines = plt.plot(xbins[1:], mu, zorder=1, label=xx[0])
#        plt.fill_between(xbins[1:], mu - std, mu + std, facecolor=lines[0].get_c(), zorder=2, alpha=0.5, lw=1)
    plt.xlabel('Average difference between $E[\\Psi]$ and group median')
    plt.ylabel('Cumulative density')
    #plt.xlim((0,0.35))
    #plt.ylim((0.6,1))
    plt.legend()
    plt.savefig(f'{args.o}', transparent=True)


if __name__ == '__main__':
    main()
        
