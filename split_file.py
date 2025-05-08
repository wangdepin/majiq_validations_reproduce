import pandas
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('filename', type=str)
parser.add_argument('-c','--columns', nargs='+', dest='columns')
parser.add_argument('-o','--output', dest='output')

args = parser.parse_args()

csv = pandas.read_csv(args.filename, sep='\t', header=0, skipinitialspace=True)
# header=3 because your header is on the third line
csv_vals = csv[args.columns]
csv_vals.to_csv(args.output, index=True, sep='\t', na_rep='nan')
