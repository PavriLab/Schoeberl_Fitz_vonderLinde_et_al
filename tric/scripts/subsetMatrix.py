import argparse as ap
import numpy as np
import logging

def get_bin_index(siteStart, leftBound, rightBound, binsize):
    binbounds = np.arange(leftBound, rightBound, binsize)
    # -1 because 0-based indices
    return len(np.where(binbounds < siteStart)[0]) - 1 \
           if not (siteStart < binbounds[0] or siteStart > binbounds[-1]) \
           else None


logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)
parser = ap.ArgumentParser()
parser.add_argument('-i', '--input')
parser.add_argument('-r', '--region')
parser.add_argument('-sr', '--subregion')
parser.add_argument('-bs', '--binsize', type = int, default = 1000)
parser.add_argument('-o', '--output')
args = parser.parse_args()

chrom1, region1 = args.region.split(':')
start1, end1 = [int(i) for i in region1.split('-')]

chrom2, region2 = args.subregion.split(':')
start2, end2 = [int(i) for i in region2.split('-')]

startbin = get_bin_index(start2, start1, end1, args.binsize)
endbin = get_bin_index(end2, start1, end1, args.binsize)

input = np.loadtxt(args.input, delimiter = '\t')
output = input[startbin:endbin + 1, startbin:endbin + 1]
np.savetxt(args.output, output, delimiter = '\t')