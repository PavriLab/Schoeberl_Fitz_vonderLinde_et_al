import argparse as ap
import matplotlib.pyplot as plt
import logging
import pandas as pd


logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)
parser = ap.ArgumentParser()
parser.add_argument('-d', '--csvs', nargs = '+',
                    help = 'space-separated list of output folder from ccseq pipeline')
parser.add_argument('-s', '--sampleNames', nargs='*',
                    help = 'optional, space-separated list of sample names corresponding to input files or named according to input directory if not specified')
parser.add_argument('--fig_width', default = 12, type = int,
                    help = 'width of the figure in inches')
parser.add_argument('-o', '--outFile', required = True,
                    help = 'file to save the plot to')
args = parser.parse_args()

twoway, threeway, nway = [], [], []
for csv in args.csvs:
    dp = pd.read_csv(csv, header=0)
    twoway.append(sum(dp['2way']) / dp.shape[0])
    threeway.append(sum(dp['3way']) / dp.shape[0])
    nway.append(sum(dp['>3way']) / dp.shape[0])
    

fig, ax = plt.subplots(1, 1)

# generate overall barchart
pos = list(range(1, len(args.sampleNames) + 1))
bottom = [0] * len(args.sampleNames)
for w, label, color in zip([nway, threeway, twoway], ['>3way', '3way', '2way'], ['darkgreen', 'limegreen', 'dimgrey']):
    ax.bar(pos, w, bottom = bottom, label = label, color = color)
    bottom = [i + j for i, j in zip(w, bottom)]

ax.set_title('n-way read fractions')
ax.legend()
ax.set_xticks(pos)
ax.set_xticklabels(args.sampleNames, rotation = 45, ha = 'right')
ax.set_ylabel('avg reads')

fig.set_figwidth(args.fig_width)
fig.set_figheight(args.fig_width)
fig.tight_layout()
fig.savefig(args.outFile)