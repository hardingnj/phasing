__author__ = 'Nicholas Harding'
__version__ = 'v1.0.0'

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


def find_reference_gaps(genome_fa):

    # fix on a gapsize of 1kb, as ld breaks down over this distance and
    # no chance of read based phasing
    # do chunks of 10k at a time. look for consecutive Ns.

    size = 10000
    gap_size = 1000

    gaps = list()

    for i in range(0, len(genome_fa), size - (2 * gap_size)):
        text = genome_fa[i:i+size]

        for m in re.finditer("N{1000}?", text, flags=re.IGNORECASE):
            gaps.append((i + m.start(), i + m.end()))

    tree = IntervalTree.from_tuples(gaps)
    tree.merge_overlaps()
    tree = sorted(tree)

    return np.array([[iv.begin, iv.end] for iv in tree])


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


# plots the colourgrid showing the position of the errors.
# Need to add a start/stop parameter as insufficient pixels to show data.
def plot_switch_errors(pos, switch_data, start=None, stop=None):

    parent_ids = list(switch_data.keys())
    parent_ids.sort(reverse=True)
    calldata = list()

    for p in parent_ids:
        switches, pos_switches, eval_pos = switch_data[p]
        c = determine_err_locs(pos, eval_pos, pos_switches[:-1].cumsum())
        calldata.append(c)
    calldata = np.array(calldata).T

    if start is not None:
        index = loc.locate_interval(pos, start, stop)
        pos = pos[index]
        calldata = calldata[index]
    else:
        start = 1
        stop = pos[-1]

    # find ratios of 1/2s
    hets = np.sum(calldata > 0, axis=1)
    errs = np.sum(calldata == 2, axis=1)
    bins = np.linspace(start, stop, int(pos.size/1000))

    midpoints = np.array([np.mean([x, y]) for x, y in zip(bins[:-1], bins[1:])])

    esites, _, _ = binned_statistic(pos, hets, np.sum, bins=bins)
    errors, _, _ = binned_statistic(pos, errs, np.sum, bins=bins)
    err_rate = errors/esites

    err_rate[esites < 10] = np.NaN

    fig = plt.figure(figsize=(12, 10))

    # plot genotypes
    ax = plt.subplot2grid((8, 1), (0, 0), rowspan=6)
    gt.plot_discrete_calldata(calldata, labels=[s.decode() for s in parent_ids],
                              colors=('grey', 'white', 'red'),
                              states=(0, 1, 2), ax=ax)

    ax = plt.subplot2grid((8, 1), (6, 0), rowspan=1)
    loc.plot_variant_locator(pos, step=pos.size/100, ax=ax,
                             line_args={"c": "k"})
    ax.xaxis.set_ticklabels([])

    ax = plt.subplot2grid((8, 1), (7, 0), rowspan=1)
    ax.plot(midpoints, err_rate, 'r')
    ax.set_xlim((start, stop))
    ul = ceil(10 * float(np.nanmax(err_rate)))/10
    ax.set_ylim((0, ul))
    ax.set_yticks(np.arange(0, ul + 0.01, 0.1))
    ax.grid(True)
    ax.set_ylabel("Switch error rate")
    ax.set_xlabel("Genomic position (bp)")

    # ax.xaxis.set_major_formatter(
    #     ticker.FuncFormatter(lambda y, ps: '%.1f' % (y*1e-6)))
    fig.subplots_adjust(hspace=0.3)
    return fig


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


def calculate_switch_distances(windows, switch_array, positions, rohz, gaps):

    positions = allel.SortedIndex(positions)
    gap_mp = np.mean(gaps, axis=1)

    marker_count = np.zeros(windows.shape[0], dtype="int")
    marker_dist = np.zeros(windows.shape[0], dtype="float")
    error_count = np.zeros(windows.shape[0], dtype="int")

    pos_sw = ph.switch.derive_position_switch_array(switch_array)
    pos_errors = np.take(positions, pos_sw[:-1].cumsum())

    for i, (start, stop) in enumerate(windows):

        # this is the code I need to change
        # A don't count error if immediately after GAP
        # B don't count towards distance
        try:
            ix = positions.locate_range(start, stop)

        except KeyError:
            marker_dist[i] = 0.0
            marker_count[i] = 0
            error_count[i] = 0
            continue

        # how many separate gaps between first and last ix?
        gap_ix = np.searchsorted(positions[ix], gap_mp)

        # interested in number of gaps
        gap_pos = np.unique(
            np.compress((gap_ix < positions[ix].size) & (gap_ix > 0), gap_ix))

        # now insert 0 and pos size at beginning and end
        cuts = np.concatenate([[0], gap_pos, [positions[ix].size]])
        assert cuts.size >= 2

        for p, q in zip(cuts[:-1], cuts[1:]):
            error_count[i] += np.sum(evaluate_markers(positions[ix][p:q],
                                                      pos_errors))

            marker_dist[i] += calc_marker_dist(positions[ix][p:q], rohz)
            # just one marker is not informative.
            marker_count[i] += (q - p - 1)

    return np.vstack([marker_dist, marker_count, error_count])


def draw_msd_across_genome(genome_pos, marker_dist, mean_switch_dist,
                           lower, upper, number_markers):

    mean_switch_dist[number_markers < 50] = np.NaN
    fig = plt.figure(figsize=(12, 3))
    ax = fig.add_subplot(111)

    ax.plot(genome_pos, 2 * marker_dist, 'k-', linewidth=2.0)
    ax.plot(genome_pos, mean_switch_dist, 'r-',  linewidth=2.0)

    # cords = np.array([[genome_pos, genome_pos[::-1],
    #                    upper, lower[::-1]]]).T.reshape((-1, 2), order="A")
    #
    # poly = plt.Polygon(cords, facecolor="r", edgecolor=None, alpha=0.2)
    #
    # ax.add_patch(poly)

    ax.grid(True)
    ax.set_ylabel("Distance (bp)")
    ax.set_xlabel("Genomic position (bp)")
    ax.set_yscale("log")
    ax.set_ylim((100, 100000))

    return fig


def draw_switch_err_rate(genome_pos, rate, rate_l, rate_u, number_markers):

    rate[number_markers < 50] = np.NaN
    rate_l[number_markers < 50] = np.NaN
    rate_u[number_markers < 50] = np.NaN

    fig = plt.figure(figsize=(12, 3))
    ax = fig.add_subplot(111)

    ax.plot(genome_pos, 1 - rate, 'r-', linewidth=2.0)
    #ax.plot(genome_pos, mean_switch_dist, 'r-',  linewidth=2.0)

    # upper = 1 - rate_u
    # lower = 1 - rate_l
    #
    # cords = np.array([[genome_pos, genome_pos[::-1],
    #                    upper, lower[::-1]]]).T.reshape((-1, 2), order="A")
    # poly = plt.Polygon(cords, facecolor="r", edgecolor=None, alpha=0.2)
    # ax.add_patch(poly)

    ax.grid(True)
    ax.set_ylim((0.82, 1.0))
    ax.set_ylabel("1 - (switch error rate)")
    ax.set_xlabel("Genomic position (bp)")

    return fig


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
t_gt = t_gt.compress(keep_pos, axis=0)

# conversely all eval MUST be in test set.
assert e_gt.shape == t_gt.shape, ("Not same shape:", e_gt.shape, t_gt.shape)

is_missing = e_gt.is_missing()

scan_windows = allel.stats.window.position_windows(pos=None, start=1,
                                                   stop=contig_length,
                                                   size=args.winsize,
                                                   step=args.overlap)

res = np.empty((len(sample_names), 3, scan_windows.shape[0]))

for idx, sid in enumerate(sample_names):

    # only consider hets in the evaluation set. ie exclude missing + homs
    hz = gt.is_het(e_gt[:, idx])
    sample_pos = np.compress(hz, e_pos)
    sample_gt = np.compress(hz, t_gt[:, idx], axis=0)
    sample_gs = np.compress(hz, e_gt[:, idx], axis=0)

    switch = ph.switch.determine_haplotype_switch(sample_gs, sample_gt)

    res[idx, :, :] = calculate_switch_distances(windows=scan_windows,
                                                switch_array=switch,
                                                positions=sample_pos,
                                                rohz=roh[sid],
                                                gaps=reference_gaps)

# now summarize across the samples dimension.
sum_r = res.sum(axis=0)

distance, n_markers, n_errors = sum_r
n_markers = n_markers.astype("int")
n_errors = n_errors.astype("int")
p_error = beta.ppf(0.5, 1 + n_errors, 1 + n_markers - n_errors)
p_error_l = beta.ppf(0.025, 1 + n_errors, 1 + n_markers - n_errors)
p_error_u = beta.ppf(0.975, 1 + n_errors, 1 + n_markers - n_errors)
err_rate = n_errors/(n_markers + n_errors)

mean_marker_d = distance/n_markers
mean_switch_d = mean_marker_d/p_error
mean_switch_d_l = mean_marker_d/p_error_l
mean_switch_d_u = mean_marker_d/p_error_u

df = pd.DataFrame.from_items((("start", scan_windows.T[0]),
                              ("stop", scan_windows.T[1]),
                              ("n_markers", n_markers),
                              ("n_errors", n_errors),
                              ("err_rate", err_rate),
                              ("distance", distance),
                              ("mean_marker_dist", mean_marker_d),
                              ("mean_switch_dist", mean_switch_d),
                              ("mean_switch_dist_2.5%", mean_switch_d_l),
                              ("mean_switch_dist_97.5%", mean_switch_d_u)))

het_positions = gt.is_het(e_gt).any(axis=1)

# 29/2/16 this plot is deprecated currently. Does not make main paper.
#pl = plot_switch_errors(np.compress(het_positions, e_pos), data)
#pl.savefig(filestem + "switch_errors.png", bbox_inches="tight")

pl = draw_msd_across_genome(scan_windows.mean(1), mean_marker_d, mean_switch_d,
                            mean_switch_d_l, mean_switch_d_u, n_markers)
pl.savefig(filestem + "mean_switch_distance.png", bbox_inches="tight")

pl = draw_switch_err_rate(scan_windows.mean(1), err_rate, p_error_l,
                          p_error_u, n_markers)
pl.savefig(filestem + "switch_error_rate.png", bbox_inches="tight")

df.to_csv(filestem + "table_switch_errors.txt", sep="\t", index=False)