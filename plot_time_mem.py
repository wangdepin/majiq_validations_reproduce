import argparse
import matplotlib.pyplot as plt
from tools import operator, all_stats
import pickle
import os

__author__    = "Jorge D. Vaquero"
__copyright__ = "Licensed under GNU GPLv3"


def plot_time(time_fname, suffix, title):

    tools = {}
    plt.figure()
    stamps = set()
    with open(time_fname) as fp:
        for ln in fp.readlines():
            tab = ln.strip().split()
            print(tab)
            try:
                tools[tab[0]][tab[1]] = get_sec(tab[2])
            except KeyError:
                tools[tab[0]] = {}
                tools[tab[0]][tab[1]] = get_sec(tab[2])
            stamps.add(tab[1])

    stamps = sorted(list(stamps))

    for kk, vv in tools.items():
        y = [vv[x] for x in stamps]

        plt.plot([1,3,6], y, label=kk)


    plt.legend()
    plt.title("Elapsed Time %s" % title)
    plt.savefig('./time_%s.pdf' % suffix)
    # plt.show()

def plot_mem(mem_fname, suffix, title):

    tools = {}
    plt.figure()
    stamps = set()
    with open(mem_fname) as fp:
        for ln in fp.readlines():
            tab = ln.strip().split()
            val = float(tab[2])/float(10**6)
            print(val, tab[2])
            try:
                tools[tab[0]][tab[1]] = val
            except KeyError:
                tools[tab[0]] = {}
                tools[tab[0]][tab[1]] = val
            stamps.add(tab[1])

    stamps = sorted(list(stamps))

    for kk, vv in tools.items():
        y = [vv[x] for x in stamps]
        print(y)
        plt.plot([1,3,6], y, label=kk)


    plt.legend()
    plt.title("Memory peak %s" % title)
    plt.savefig('./mem_%s.pdf' % suffix)
    # plt.show()


def get_sec(time_str):
    """Get Seconds from time."""
    h, m, s = time_str.split(':')
    return int(h) * 3600 + int(m) * 60 + int(s)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('time_input')
    parser.add_argument('mem_input')
    parser.add_argument('--suffix', action='store', dest='suffix' )
    parser.add_argument('--title', action='store', dest='title' )

    args = parser.parse_args()

    plot_time(args.time_input, args.suffix, args.title)
    plot_mem(args.mem_input, args.suffix, args.title)
