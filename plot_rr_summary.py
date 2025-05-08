import matplotlib
matplotlib.use('agg')
import argparse
from tools import operator, all_stats
import pickle
import os
import numpy as np
import matplotlib.pyplot as plt
import errno


def main(input):

    data_rr = {}
    data_n = {}
    x_axis = set()
    with open(input) as fp:
        for li in fp.readlines():
            if li.startswith('#'):
                continue
            tab = li.strip().split()
            x_id = int(tab[1].split('vs')[0])
            try:
                data_n[tab[0]][x_id] = int(tab[2])
                data_rr[tab[0]][x_id] = float(tab[3])
            except KeyError:
                data_n[tab[0]] = {x_id: int(tab[2])}
                data_rr[tab[0]] = {x_id: float(tab[3])}
            x_axis.add(x_id)


    plt.figure(figsize=(16, 4))
    x_axis = sorted(list(x_axis))
    for xx, vv in data_rr.items():
        X = []
        Y = []
        for y in x_axis:
            if y in vv:
                X.append(y)
                Y.append(vv[y])

        plt.plot(X, Y, label=xx, marker='o')

    plt.legend()
    plt.savefig('summary_rr.pdf')

    plt.figure(figsize=(16, 4))
    x_axis = sorted(list(x_axis))
    for xx, vv in data_n.items():
        X = []
        Y = []
        for y in x_axis:
            if y in vv:
                X.append(y)
                Y.append(vv[y])

        plt.plot(X, Y, label=xx, marker='o')

    plt.legend()
    plt.savefig('summary_rrN.pdf')



if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('input')
    args = parser.parse_args()

    main(args.input)