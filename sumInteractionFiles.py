import argparse as ap
import pandas as pd
import logging

logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)
parser = ap.ArgumentParser()
parser.add_argument('-i', '--input', nargs = '+',
                    help = 'probe interaction files')
parser.add_argument('-o', '--outputFile', required = True,
                    help = 'name of the output file')
args = parser.parse_args()

tab = pd.DataFrame(columns = ['coords1', 'coords2', 'count'])
for filename in args.input:
    tmp_tab = pd.read_csv(filename,
                          sep = '\t',
                          header = None,
                          names = ['coords1', 'coords2', 'count'])
    tab = tab.merge(tmp_tab, on = ['coords1', 'coords2'], how = 'outer') \
             .fillna(0)

tab['count'] = tab.filter(regex = 'count_*') \
                  .sum(axis = 1)

print(tab['count'].max())
tab[['coords1', 'coords2', 'count']].to_csv(args.outputFile,
                                            sep = '\t',
                                            header = False,
                                            index = False)