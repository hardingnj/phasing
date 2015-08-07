#! /usr/bin/python
from __future__ import print_function

__author__ = 'Nicholas Harding'

import phasing as ph
import h5py
import numpy as np
import anhima
import argparse
import time
import hmmlearn
from hmmlearn import hmm
import re

check = re.compile(".+\{chrom\}.+")

start_time = time.time()
chunk_size = 2e5
parser = argparse.ArgumentParser(
    description='Tool to identify runs of homozygosity from sequence data')

# files:
parser.add_argument('input', help='input hdf5 filestem')
parser.add_argument('output', help='output file')
parser.add_argument('--accessibility', '-A', action='store', default=None,
                    dest='accessibility', help='Accessibility h5',
                    required=True)

# data:
parser.add_argument('--sample', '-S', action='store', default=None,
                    dest='sample', required=True,
                    help='Which sample to evaluate.')
parser.add_argument('--chr', '-C', nargs='+',
                    default=["2L", "2R", "3L", "3R", "X"],
                    dest='chrom', help='Which contig to model')

# parameters
parser.add_argument('--emission', '-X', action='store',
                    default=[0.0001, 0.001, 0.015],
                    dest='phet', required=True, type=float, nargs="+",
                    help='Probability of a het in a run of homozygozity.')

parser.add_argument('--transition', '-T', action='store', default=1e-6,
                    type=float, dest='transition',
                    help='Transition probability of model')

args = parser.parse_args()

assert check.match(args.input), \
    "Input file does not have format string of {chrom}"


def predict_roh_state(model, ind_genotype, pos, accessible, label="unk"):

    # assume non-included pos are homozygous reference
    assert ind_genotype.shape[0] == pos.size

    heterozygosity = anhima.gt.is_het(ind_genotype)
    observations = np.zeros(accessible.size, dtype='int')
    for i in np.compress(heterozygosity, pos):
        observations[i - 1] = 1

    for i in np.where(np.invert(accessible))[0]:
        observations[i - 1] = 2

    predictions = model.predict(obs=observations)
    print('Predictions complete for {0}:'.format(label))
    print('{0}% ROH'.format(100*np.mean(predictions == 0)))
    print('{0}% Normal'.format(100*np.mean(predictions >= 1)))

    return predictions, observations


# This function is slightly hokey. Relies on there being > 1 of a type.
def get_state_windows(predicted_state, state=0):

    assert isinstance(predicted_state, np.ndarray), \
        "get_state_windows expects an ndarray"
    wh = np.where(predicted_state == state)[0]
    if wh.size == 0:
        # then there are no things of interest
        return None
    elif wh.size == 1:
        # the whole thing is one big thing of interest
        return list(np.array([1, predicted_state.size]))

    intervals = list()
    iv_start = wh[0]

    pos = None
    for i, pos in enumerate(wh[1:]):
        if (pos - wh[i]) > 1:
            intervals.append([iv_start, wh[i]])
            iv_start = pos
    intervals.append([iv_start, pos])

    roh = np.array(intervals)
    print("{0} distinct windows of state={1}".format(roh.shape[0], state))
    # correct for fact that pos 1 is in index 0.
    return roh + 1


def derive_emission_mx(prob, pa):
    # one row per p in prob
    # hom, het, unobserved
    mx = [[(1-p) * (1-pa), p * (1-pa), pa] for p in prob]
    mxe = np.array(mx)
    assert mxe.shape == (prob.size, 3)
    return mxe


def derive_transition_mx(pr, nstates):
    # this is a symmetric matrix
    mx = np.zeros((nstates, nstates))
    effective_tp = pr/(nstates-1)
    for i in range(nstates):
        for j in range(nstates):
            if i == j:
                mx[i, j] = 1 - pr
            else:
                mx[i, j] = effective_tp
    return mx


def calculate_windows(contig):
    # this is the main function: given a chrom and parameters, it computes the
    # likely ROH using an HMM model.

    snps_fn = args.input.format(chrom=contig)
    print(snps_fn)
    fh = h5py.File(snps_fn, "r")[contig]
    acc = h5py.File(args.accessibility, "r")[contig]

    is_accessible = acc['is_accessible'][:]

    # data
    samples = fh["samples"][:]
    positions = fh['variants']['POS'][:]

    # define model:
    ### 2 state model, 3 observable
    # 0: ROH
    # 1: normal
    # 2: inaccessible
    emission_px = np.array(args.phet)
    emission_px.sort()

    # start probability
    start_prob = np.repeat(1/emission_px.size, emission_px.size)

    # transition probabilty:
    trans_p = args.transition

    # transition between underlying states
    transitions = derive_transition_mx(trans_p, emission_px.size)

    # probability of inaccessible
    p_accessible = is_accessible.mean()

    emission_mx = derive_emission_mx(emission_px, p_accessible)

    # initialize HMM
    roh_hmm = hmm.MultinomialHMM(n_components=emission_px.size)

    roh_hmm.n_symbols_ = 3
    roh_hmm.startprob_ = start_prob
    roh_hmm.transmat_ = transitions
    roh_hmm.emissionprob_ = emission_mx

    idx = samples.tolist().index(args.sample.encode())
    genotype = fh['calldata']['genotype'][:, idx]

    pred, obs = predict_roh_state(roh_hmm, genotype, positions,
                                  is_accessible, label=contig)

    #inaccessible_windows = get_state_windows(obs, state=2)
    homozygous_windows = get_state_windows(pred, state=0)

    return homozygous_windows

# dump settings to file.
with open(args.output, "w") as out:
    for k, v in vars(args).items():
        print("# {0}: {1}".format(k, v), file=out)
    print("# v. phasing: " + ph.__version__, file=out)
    print("# v. h5py: " + h5py.__version__, file=out)
    print("# v. numpy: " + np.__version__, file=out)
    print("# v. anhima: " + anhima.__version__, file=out)
    print("# v. argparse: " + argparse.__version__, file=out)
    print("# v. hmmlearn: " + hmmlearn.__version__, file=out)
    print("# v. re: " + re.__version__, file=out)

for chrom in args.chrom:
    homozygous = calculate_windows(chrom)
    with open(args.output, "a") as out:
        for w in homozygous:
            print("\t".join([chrom] + w.astype('str').tolist()), file=out)

print("---------------------------")
print("Completed in {0} minutes".format((time.time() - start_time)/60))