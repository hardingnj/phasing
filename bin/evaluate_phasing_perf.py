__author__ = 'Nicholas Harding'
__version__ = 'v1.2.0'

# This script takes a list of samples,
# two sets of phased data and evaluates one against the other

# to produce:
# - 2 plots of performance.
    # - A raw data file describing all switch errors
# - A summary table of number of errors in each window.

import argparse
import h5py
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import allel
import anhima.gt as gt
import anhima.loc as loc
import phasing as ph
from scipy.stats import binned_statistic, beta
from os.path import join
from math import ceil
import pyfasta
import re
from intervaltree import IntervalTree


def haldane(r):
    return -np.log(1-(2*r))/2


def find_reference_gaps(genome_fa):

    # fix on a gapsize of 1kb, as ld breaks down over this distance and
    # no chance of read based phasing
    # do chunks of 10k at a time. look for consecutive Ns.
    size = 10000
    gap_size = 1000
    match = "N{{{n}}}?".format(n=gap_size)

    gaps = list()

    for i in range(0, len(genome_fa), size - (2 * gap_size)):
        text = genome_fa[i:i+size]

        for m in re.finditer(match, text, flags=re.IGNORECASE):
            gaps.append((i + m.start(), i + m.end()))

    tree = IntervalTree.from_tuples(gaps)
    tree.merge_overlaps()

    return np.array([[iv.begin, iv.end] for iv in sorted(tree)])


# utility function for plot_switch_errors,
# classifies each position as het good, het bad or non informative.
def determine_err_locs(positions, evaluated_positions, switches_indexes ):

    position_switches = np.take(evaluated_positions, switches_indexes)
    result = np.zeros(positions.shape)

    for pos in evaluated_positions:
        r = 1
        if pos in position_switches:
            r = 2
        result[loc.locate_position(positions, pos)] = r
    return result


def calc_marker_dist(pos, homozygosity=None):

    if pos.size == 0:
        return 0
    else:
        start, stop = pos[0], pos[-1]

    if homozygosity is None:
        sum_roh = 0
    else:
        int_tree = IntervalTree(IntervalTree.from_tuples(homozygosity))
        int_tree.slice(start)
        int_tree.slice(stop)
        sum_roh = np.sum([ival.length() for ival in int_tree[start:stop]])

    return stop - start - sum_roh


def evaluate_markers(markers, error_positions):
    """
    function that retuns all marker distances, and whether they are a SE.
    start: is the genomic window start pos
    markers: an arr of genetic positions of the markers
    error_pos: an arr of error positions
    window_size: window size
    """

    if markers.size == 0:
        return np.empty(0, dtype="bool")

    # don't care if first marker is an error as window based ie pertains to
    # gap between marker and immediately before
    errors = np.in1d(markers[1:], error_positions)

    return errors


def calculate_switch_distances(windows, switch_array, marker_pos, hz_pos,
                               rohz, gaps):

    marker_pos = allel.SortedIndex(marker_pos)
    gap_mp = np.mean(gaps, axis=1)

    assert np.in1d(marker_pos, hz_pos).all(), "all markers are subset of hets"

    marker_count = np.zeros(windows.shape[0], dtype="int")
    marker_dist = np.zeros(windows.shape[0], dtype="float")
    error_count = np.zeros(windows.shape[0], dtype="int")
    hz_count = np.zeros(windows.shape[0], dtype="int")

    pos_sw = ph.switch.derive_position_switch_array(switch_array)
    pos_errors = np.take(marker_pos, pos_sw[:-1].cumsum())

    for i, (start, stop) in enumerate(windows):

        # this is the code I need to change
        # A don't count error if immediately after GAP
        # B don't count towards distance
        try:
            ix = marker_pos.locate_range(start, stop)

        except KeyError:
            marker_dist[i] = 0.0
            marker_count[i] = 0
            error_count[i] = 0
            hz_count[i] = 0
            continue

        # how many separate gaps between first and last ix?
        gap_ix = np.searchsorted(marker_pos[ix], gap_mp)

        # interested in number of gaps
        gap_pos = np.unique(
            np.compress((gap_ix < marker_pos[ix].size) & (gap_ix > 0), gap_ix))

        # now insert 0 and pos size at beginning and end
        cuts = np.concatenate([[0], gap_pos, [marker_pos[ix].size]])
        assert cuts.size >= 2

        for p, q in zip(cuts[:-1], cuts[1:]):

            first, last = marker_pos[ix][p], marker_pos[ix][q-1]

            # how many hets between first an last?
            counthets = np.searchsorted(hz_pos, last) - \
                np.searchsorted(hz_pos, first)

            error_count[i] += np.sum(evaluate_markers(marker_pos[ix][p:q],
                                                      pos_errors))

            marker_dist[i] += calc_marker_dist(marker_pos[ix][p:q], rohz)
            # just one marker is not informative.
            marker_count[i] += (q - p - 1)
            hz_count[i] += counthets

    return np.vstack([marker_dist, marker_count, error_count, hz_count])


parser = argparse.ArgumentParser(
    description='Evaluation pipeline for phasing evaluation')

# data:
parser.add_argument('--test', '-T', help='input hdf5 file', action='store',
                    dest="test", default=None, type=str)

parser.add_argument('--eval', '-E', help='input hdf5 file to evaluate against',
                    action='store', dest="eval", default=None, type=str)

parser.add_argument('--output', '-O', help='output directory', action="store",
                    dest="outdir", default="/tmp/", type=str)

parser.add_argument('--fasta', '-F', action='store', default=None,
                    dest='fasta', help='FASTA file', type=str,
                    required=True)

parser.add_argument('--accessibility', '-A', action='store', default=None,
                    dest='accessibility', help='Accessibility h5', type=str,
                    required=False)

parser.add_argument('--samples', '-S', action='store', nargs="+", default=None,
                    dest='samples', required=False,
                    help='Which samples to evaluate.')

parser.add_argument('--chr', '-C', default=None, required=True, action="store",
                    dest='chrom', help='Which contig to evaluate')

parser.add_argument('--stem', default=None, dest='stem',
                    help='filestem')

# non-required parameters
parser.add_argument('--roh', '-R', action='store',
                    default=None, dest='roh', type=str,
                    help='Path to ROH file for calculating switch distance')

parser.add_argument('--overlap', '-wo', action='store',
                    default=250000, dest='overlap', type=int,
                    help='Overlap between windows.')

parser.add_argument('--windowsize', '-ws', action='store', default=500000,
                    type=int, dest='winsize',
                    help='Evenly spaced windows across genome')


args = parser.parse_args()
genome = pyfasta.Fasta(args.fasta)
contig_length = len(genome[args.chrom])
reference_gaps = find_reference_gaps(genome[args.chrom])

filestem = join(args.outdir, "{chrom}_{stem}_".format(chrom=args.chrom,
                                                      stem=args.stem))

test_fh = h5py.File(args.test, "r")
eval_fh = h5py.File(args.eval, "r")

eval_samples = eval_fh[args.chrom]["samples"][:].tolist()
test_samples = test_fh[args.chrom]['samples'][:].tolist()


if args.samples is not None:
    test_idx = [test_samples.index(s.encode()) for s in args.samples]
    eval_idx = [eval_samples.index(s.encode()) for s in args.samples]
    sample_names = [s.encode() for s in args.samples]
else:
    intersection = np.intersect1d(eval_samples, test_samples)
    test_idx = [test_samples.index(s) for s in intersection]
    eval_idx = [eval_samples.index(s) for s in intersection]
    sample_names = intersection

if args.roh is not None:
    roh_fh = h5py.File(args.roh, "r")
    roh = dict()
    for s in sample_names:
        try:
            roh[s] = roh_fh[args.chrom][s][:]
        except KeyError:
            roh[s] = None
else:
    roh = {s: None for s in sample_names}


t_geno = allel.GenotypeCArray.from_hdf5(test_fh[args.chrom]["calldata"][
    "genotype"])
t_gt = t_geno.take(test_idx, axis=1)
e_geno = allel.GenotypeCArray.from_hdf5(eval_fh[args.chrom]["calldata"][
    "genotype"])
e_gt = e_geno.take(eval_idx, axis=1)

# now lets line up postions...
e_pos = eval_fh[args.chrom]['variants']['POS'][:]
t_pos = test_fh[args.chrom]['variants']['POS'][:]

# keep only the test positions that are in the eval set.
keep_pos = np.in1d(t_pos, e_pos)
reduced_test_gt = t_gt.compress(keep_pos, axis=0)

# conversely all eval MUST be in test set.
assert e_gt.shape == reduced_test_gt.shape, \
    ("Not same shape:", e_gt.shape, reduced_test_gt.shape)

is_missing = e_gt.is_missing()

scan_windows = allel.stats.window.position_windows(pos=None, start=1,
                                                   stop=contig_length,
                                                   size=args.winsize,
                                                   step=args.overlap)

res = np.empty((len(sample_names), 4, scan_windows.shape[0]))

# all the hz in the test set
het_pos = t_gt.is_het()

for idx, sid in enumerate(sample_names):

    # only consider hets in the evaluation set. ie exclude missing + homs
    hz = gt.is_het(e_gt[:, idx])
    sample_pos = np.compress(hz, e_pos)
    sample_gt = np.compress(hz, reduced_test_gt[:, idx], axis=0)
    sample_gs = np.compress(hz, e_gt[:, idx], axis=0)

    switch = ph.switch.determine_haplotype_switch(sample_gs, sample_gt)

    all_het_pos = np.compress(het_pos[:, idx], t_pos, axis=0)

    res[idx, :, :] = calculate_switch_distances(windows=scan_windows,
                                                switch_array=switch,
                                                marker_pos=sample_pos,
                                                hz_pos=all_het_pos,
                                                rohz=roh[sid],
                                                gaps=reference_gaps)

# now summarize across the samples dimension.
sum_r = res.sum(axis=0)

distance, n_markers, n_errors, n_hets = sum_r
n_markers = n_markers.astype("int")
n_errors = n_errors.astype("int")
n_hets = n_hets.astype("int")

marker_err_rate = n_errors/(n_markers + n_errors)

# estimate the per-het error rate
het_err_rate = haldane(marker_err_rate)/(n_hets/n_markers)

mean_marker_d = distance/n_markers
mean_switch_d = mean_marker_d/marker_err_rate

df = pd.DataFrame.from_items((("start", scan_windows.T[0]),
                              ("stop", scan_windows.T[1]),
                              ("n_hets", n_hets),
                              ("n_markers", n_markers),
                              ("n_errors", n_errors),
                              ("marker_err_rate", marker_err_rate),
                              ("het_err_rate", het_err_rate),
                              ("distance", distance),
                              ("mean_marker_dist", mean_marker_d),
                              ("mean_switch_dist", mean_switch_d)))

het_positions = gt.is_het(e_gt).any(axis=1)

df.to_csv(filestem + "table_switch_errors.txt", sep="\t", index=False)