#! /usr/bin/python
from __future__ import print_function

__author__ = 'Nicholas Harding'

import re
import phasing as ph
import h5py
import numpy as np
import anhima
import argparse
import time
from hmmlearn import hmm

start_time = time.time()
chunk_size = 2e5
parser = argparse.ArgumentParser(
    description='Tool to identify runs of homozygosity from sequence data')

parser.add_argument('input', help='input hdf5 file')
parser.add_argument('output', help='output file')

parser.add_argument('--sample', '-S', action='store', default=None,
                    dest='sample', help='Which sample to evaluate.')
parser.add_argument('--chr', '-C', action='store', default=None,
                    dest='chrom', help='Which contig')
parser.add_argument('--accessibility', '-A', action='store', default=None,
                    dest='accessibility', help='Accessibility h5')

args = parser.parse_args()


def predict_roh_state(model, ind_genotype, pos, accessible):

    # assume non-included pos are homozygous reference
    assert ind_genotype.shape[0] == pos.size

    heterozygosity = anhima.gt.is_het(ind_genotype)
    observations = np.zeros(accessible.size, dtype='int')
    for i in np.compress(heterozygosity, pos):
        observations[i] = 1

    for i in np.where(~accessible)[0]:
        observations[i] = 2

    predictions = model.predict(obs=observations)
    print('Predictions complete:')
    print('{0}% ROH'.format(100*np.mean(predictions == 0)))
    print('{0}% Normal'.format(100*np.mean(predictions == 1)))

    return predictions, observations


def get_state_windows(predicted_state, state=0):

    x = np.where(np.diff(predicted_state == state))[0]

    if predicted_state[0] == state:
        x = np.concatenate([[0], x])
    if predicted_state[-1] == state:
        x = np.concatenate([x, [predicted_state.size - 1]])

    roh = x.reshape(-1, 2)
    print("{0} distinct windows of state={1}".format(roh.shape[0], state))
    return roh

# load files
fh = h5py.File(args.input, "r")[args.chrom]
acc = h5py.File(args.accessibility, "r")[args.chrom]

is_accessible = acc['is_accessible'][:]

# data
samples = fh["samples"][:]
positions = fh['variants']['POS'][:]

# define model:
### 2 state model, 3 observable
# 0: ROH
# 1: normal
# 2: inaccessible

# start probability
start_prob = np.array([0.5, 0.5])

# transition probabilty:
trans_p = 0.001

# transition between ROH and regular
transitions = np.array([[1 - trans_p, trans_p],
                        [trans_p, 1 - trans_p]])

# probability of oberving a het in ROH vs reg   ularly in an accessible region
p_het_roh = 1e-5
p_het_norm = 5e-3

# probability of inaccessible
p_inacess = is_accessible.mean()

confusion = np.array([[(1 - p_inacess) * (1 - p_het_roh),
                       (1 - p_inacess) * p_het_roh,
                       p_inacess],
                      [(1 - p_inacess) * (1 - p_het_norm),
                       (1 - p_inacess) * p_het_norm,
                       p_inacess]])

roh_hmm = hmm.MultinomialHMM(n_components=2)

roh_hmm.n_symbols_ = 3
roh_hmm.startprob_ = start_prob
roh_hmm.transmat_ = transitions
roh_hmm.emissionprob_ = confusion

idx = samples.tolist().index(args.sample)
genotype = fh['calldata']['genotype'][:, idx]

pred, obs = predict_roh_state(roh_hmm,
                              genotype,
                              positions,
                              is_accessible)

inaccessible_windows = get_state_windows(obs, state=2)
homozygous_windows = get_state_windows(pred, state=0)

with open(args.output, "w") as out:
    for k, v in vars(args).iteritems():
        print("# {0}: {1}".format(k, v), file=out)
    print("# VERSION: " + ph.__version__, file=out)
    print("# TRANSITION: " + re.sub("\n", ",", np.str(roh_hmm.transmat_)),
          file=out)
    print("# EMISSION: " + re.sub("\n", ",", np.str(roh_hmm.emissionprob_)),
          file=out)
    print("# STARTPROB: " + re.sub("\n", ",", np.str(roh_hmm.startprob_)),
          file=out)

    for w in homozygous_windows:
        print("\t".join(w.astype('str').tolist()), file=out)

print("---------------------------")
print("Completed in {0} minutes".format((time.time() - start_time)/60))