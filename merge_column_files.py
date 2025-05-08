import argparse
import os
import sys


def main(file_list, outDir, cols, cmd=""):

    if not os.path.exists(outDir):
        os.makedirs(outDir)

    merged_dict = {}

    ncols = len(file_list)
    print (ncols)
    for eidx, fle in enumerate(file_list):
        with open(fle, 'r') as fp:
            for ll in fp.readlines():
                if ll.startswith('gene'):
                    continue
                tab = ll.strip().split()
                try:
                    merged_dict[tab[0]][eidx] = tab[1]
                except KeyError:
                    merged_dict[tab[0]] = ['0'] * ncols
                    merged_dict[tab[0]][eidx] = tab[1]

    with open('%s/merged_vals.tsv' % outDir, 'w+') as fp:
        fp.write("#Command: %s\n" % cmd)
        fp.write("#ID\t%s\n" % '\t'.join(cols))
        for kk, vv in merged_dict.items():
            vals = "\t".join(vv)
            fp.write("%s\t%s\n" % (kk, vals))


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('file_list', nargs='+', help='List of files to merge, based on ID. The IDs should be equivalent:'
                                                     '\tID\tVALUE')

    parser.add_argument('-o', '--output', action='store', dest='outDir', help='Output directory where the output '
                                                                              'files will be stored')

    parser.add_argument('--cols', nargs='+', action='store', dest='cols', default=None,
                        help='Columns name to be used instead filenames')

    args = parser.parse_args()

    cmd = " ".join(sys.argv)

    if args.cols is not None:
        assert len(args.cols) == len(args.file_list), "column list and file list should have the same length"
        cols = args.cols
    else:
        cols = args.file_list



    main(args.file_list, args.outDir, cols=cols, cmd=cmd)