__author__ = 'Nicholas Harding'

# This script takes a list of samples,
# two sets of phased data and evaluates one against the other

# to produce:
# - 3 plots of performance.
# - A raw data file describing all switch errors
# - A summary table of number of errors in each window.

import argparse
import h5py
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from intervaltree import IntervalTree
import allel
import anhima.gt as gt
import anhima.loc as loc
import phasing as ph
from scipy.stats import binned_statistic, beta
from os.path import join
from math import ceil


contig_lengths = {"2L": 49364325, "2R": 61545105,
                  "3L": 41963435, "3R": 53200684,
                  "X": 24393108, "Y_unplaced": 237045,
                  "UNKN": 42389979}


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


def evaluate_markers(start, markers, error_positions, window_size=2e5):
    """
    function that retuns all marker distances, and whether they are a SE.
    start: is the genomic window start pos
    markers: an arr of genetic positions of the markers
    error_pos: an arr of error positions
    window_size: window size

    """
    sliced = loc.locate_interval(markers, start, start + window_size)
    window_markers = markers[sliced]
    if window_markers.size == 0:
        return np.empty(0, dtype="int"), np.empty(0, dtype="bool")

    # don't care if first marker is an error...
    errors = np.in1d(window_markers[1:], error_positions)
    return window_markers, errors


def calculate_switch_distances(start_points, switch_data,
                               window_size, roh_dict):

    # intial declatations
    mean_marker_dist = np.repeat(np.NaN, start_points.size)
    mean_switch_dist = np.repeat(np.NaN, start_points.size)
    conf_switch_dist = np.repeat(np.NaN, 2*start_points.size).reshape((
        start_points.size, 2))

    marker_count = np.zeros(start_points.size, dtype="int")
    error_count = np.zeros(start_points.size, dtype="int")
    pids = switch_data.keys()

    evaluated_pos = {pid: switch_data[pid][2] for pid in pids}

    # pos errors represents the position of switch errors in the chromosome.
    pos_errors = {pid: np.take(evaluated_pos[pid],
                               switch_data[pid][1][:-1].cumsum())
                  for pid in pids}

    for i, start in enumerate(start_points):

        d = {pid: evaluate_markers(start,
                                   switch_data[pid][2],
                                   pos_errors[pid],
                                   window_size) for pid in pids}

        positions = {pid: d[pid][0] for pid in pids}
        errors = {pid: d[pid][1] for pid in pids}

        distances = [calc_marker_dist(positions[pid],
                                      roh_dict[pid]) for pid in pids]

        total_distance = np.sum(distances)

        n_markers = int(np.sum([np.size(p) - 1 for p in positions.values()]))
        n_error = int(np.sum([np.sum(e) for e in errors.values()]))

        max_p_err = beta.ppf(0.975, 1 + n_error, 1 + n_markers - n_error)
        min_p_err = beta.ppf(0.025, 1 + n_error, 1 + n_markers - n_error)
        p_error = beta.ppf(0.5, 1 + n_error, 1 + n_markers - n_error)

        mean_marker_d = total_distance/n_markers
        mean_switch_d = mean_marker_d/p_error
        mean_marker_dist[i] = mean_marker_d
        mean_switch_dist[i] = mean_switch_d
        conf_switch_dist[i] = (mean_marker_d/min_p_err, mean_marker_d/max_p_err)
        marker_count[i] = n_markers
        error_count[i] = n_error

    return (mean_marker_dist, mean_switch_dist, conf_switch_dist,
            marker_count, error_count)


def draw_msd_across_genome(genome_pos, marker_dist, mean_switch_dist,
                           conf_mean_switch_dist, number_markers):

    mean_switch_dist[number_markers < 50] = np.NaN
    fig = plt.figure(figsize=(12, 3))
    ax = fig.add_subplot(111)

    ax.plot(genome_pos, 2 * marker_dist, 'k-', linewidth=2.0)
    ax.plot(genome_pos, mean_switch_dist, 'r-',  linewidth=2.0)

    upper = conf_mean_switch_dist.T[0]
    lower = conf_mean_switch_dist.T[1]

    cords = np.array([[genome_pos, genome_pos[::-1],
                       upper, lower[::-1]]]).T.reshape((-1, 2), order="A")
    poly = plt.Polygon(cords, facecolor="r", edgecolor=None, alpha=0.2)

    ax.add_patch(poly)
    ax.grid(True)
    ax.set_ylabel("Distance (bp)")
    ax.set_xlabel("Genomic position (bp)")
    ax.set_yscale("log")

    return fig


def draw_switch_err_rate(genome_pos, marker_dist, mean_switch_dist,
                         conf_mean_switch_dist, number_markers):

    mean_switch_dist[number_markers < 50] = np.NaN

    fig = plt.figure(figsize=(12, 3))
    ax = fig.add_subplot(111)

    ax.plot(genome_pos, 1 - marker_dist/mean_switch_dist, 'r-', linewidth=2.0)
    #ax.plot(genome_pos, mean_switch_dist, 'r-',  linewidth=2.0)

    upper = 1 - marker_dist/conf_mean_switch_dist.T[0]
    lower = 1 - marker_dist/conf_mean_switch_dist.T[1]

    cords = np.array([[genome_pos, genome_pos[::-1],
                       upper, lower[::-1]]]).T.reshape((-1, 2), order="A")
    poly = plt.Polygon(cords, facecolor="r", edgecolor=None, alpha=0.2)
    ax.add_patch(poly)

    ax.grid(True)
    ax.set_ylabel("1 - (switch error rate)")
    ax.set_xlabel("Genomic position (bp)")

    return fig


def switch_distances_tofile(path, mean_marker_d, mean_switch_d, conf_int, nmark,
                            nerr, window_starts, window_size):

    columns = "window", "mean marker dist", "mean switch dist", "CI 95%", \
              "n Errors", "n Markers"

    swindows = ["{0}-{1} Mb".format(np.round(x*1e-6, 2),
                                    np.round((y-1)*1e-6, 2))
                for x, y in zip(window_starts, window_starts + window_size)]

    ci = ["{0},{1}".format(np.round(x, 2), np.round(y, 2))
          for y, x in conf_int]

    frame = pd.DataFrame(columns=columns)
    frame[columns[0]] = swindows
    frame[columns[1]] = np.round(mean_marker_d, 2)
    frame[columns[2]] = np.round(mean_switch_d, 2)
    frame[columns[3]] = ci
    frame[columns[4]] = nerr
    frame[columns[5]] = nmark
    frame["error rate"] = np.round(nerr/nmark, 4)
    frame.to_csv(path + ".txt", sep="\t", index=False)
    frame.to_latex(path + ".tex", index=False)


parser = argparse.ArgumentParser(
    description='Evaluation pipeline for phasing evaluation')

# data:
parser.add_argument('--test', '-T', help='input hdf5 file', action='store',
                    dest="test", default=None, type=str)

parser.add_argument('--eval', '-E', help='input hdf5 file to evaluate against',
                    action='store', dest="eval", default=None, type=str)

parser.add_argument('--output', '-O', help='output directory', action="store",
                    dest="outdir", default="/tmp/", type=str)

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
                    default=500000, dest='overlap', type=int,
                    help='Overlap between windows.')

parser.add_argument('--windowsize', '-ws', action='store', default=250000,
                    type=int, dest='winsize',
                    help='Evenly spaced windows across genome')

args = parser.parse_args()
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

data = dict()

for idx, sid in enumerate(sample_names):

    # only consider hets in the evaluation set. ie exclude missing + homs
    hz = gt.is_het(e_gt[:, idx])
    sample_pos = np.compress(hz, e_pos)
    sample_gt = np.compress(hz, t_gt[:, idx], axis=0)
    sample_gs = np.compress(hz, e_gt[:, idx], axis=0)

    switch = ph.switch.determine_haplotype_switch(sample_gs, sample_gt)
    pos_sw = ph.switch.derive_position_switch_array(switch)
    data[sid] = switch, pos_sw, sample_pos

het_positions = gt.is_het(e_gt).any(axis=1)
pl = plot_switch_errors(np.compress(het_positions, e_pos), data)
pl.savefig(filestem + "switch_errors.png", bbox_inches="tight")

# now build the windowed MSD array.
windows = np.arange(1, contig_lengths[args.chrom], args.winsize).astype('int')

mmd, msd, conf, nmarker, nerror = calculate_switch_distances(
    windows, data, window_size=args.overlap, roh_dict=roh)

pl = draw_msd_across_genome(windows, mmd, msd, conf, nmarker)
pl.savefig(filestem + "mean_switch_distance.png", bbox_inches="tight")

pl = draw_switch_err_rate(windows, mmd, msd, conf, nmarker)
pl.savefig(filestem + "switch_error_rate.png", bbox_inches="tight")

switch_distances_tofile(filestem + "table_switch_errors", mmd, msd, conf,
                        nmarker, nerror, windows, args.overlap)