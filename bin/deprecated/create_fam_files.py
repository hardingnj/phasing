__author__ = 'Nicholas Harding'

import os
import pandas as pd
import re
import argparse

parser = argparse.ArgumentParser(description='Tool to filter a hdf5 file')
parser.add_argument('input', help='input sample pedfile')
parser.add_argument('allsamples', help='File containing all samples')
parser.add_argument('output', help='output directory')
args = parser.parse_args()

# create tfam file: contains ID, pedigree, mother father (optional sex
# and phenotype) for all crosses combined, and for individual crosses
if not os.path.isdir(args.output):
    os.mkdir(args.output)

assert os.path.isfile(args.input) and os.path.isfile(args.allsamples)

ped_tbl = pd.read_csv(args.input, sep="\t", index_col=1)
all_tbl = pd.read_csv(args.allsamples, sep="\t", header=None)
all_tbl.columns = ['Individual_ID']
all_tbl.Individual_ID = [re.sub('_', '-', s) for s in all_tbl.Individual_ID]
all_tbl['pedigree'] = pd.Series()
all_tbl['Paternal_ID'] = -1
all_tbl['Maternal_ID'] = -1
all_tbl['Sex'] = 0
all_tbl['Phenotype'] = 0

for i, row in all_tbl.iterrows():
    try:
        match = ped_tbl.loc[row['Individual_ID']]
        cross = match['cross']
        all_tbl['pedigree'][i] = cross

        parents = ped_tbl[(ped_tbl.cross == cross) & (ped_tbl.role == 'parent')].index
        if match.name in parents:
            all_tbl['Paternal_ID'][i] = '0'
            all_tbl['Maternal_ID'][i] = '0'
        else:
            all_tbl['Paternal_ID'][i] = parents[0]
            all_tbl['Maternal_ID'][i] = parents[1]
    except KeyError:
        all_tbl['pedigree'][i] = 'none'
        all_tbl['Paternal_ID'][i] = '0'
        all_tbl['Maternal_ID'][i] = '0'

# Assign the family id. Basically {family}_{counter}
di = dict()
l = list()

for _, p in all_tbl.pedigree.iteritems():
    di[p] = di.get(p, 0) + 1
    l.append(p + '_' + str(di[p]))

all_tbl['Family_ID'] = l
all_tbl[['Family_ID', 'Individual_ID', 'Paternal_ID', 'Maternal_ID', 'Sex',
         'Phenotype']].to_csv(os.path.join(args.output, 'all_crosses.fam'),
                              index=False, header=False, sep=" ")
print 'family file written'

crosses = set(ped_tbl.cross)
for c in crosses:
    ct = all_tbl[all_tbl.pedigree == c]
    ct[['Family_ID', 'Individual_ID', 'Paternal_ID', 'Maternal_ID', 'Sex',
        'Phenotype']].to_csv(os.path.join(args.output, c + '.fam'),
                             index=False, header=False, sep=" ")
    print 'Written ' + c
