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
import matplotlib.ticker as ticker
import yaml
import pandas as pd
import allel
import bcolz
import anhima.gt as gt
import anhima.loc as loc
import phasing as ph
from scipy.stats import binned_statistic, beta
from os.path import join
from math import ceil


# utility function for plot_switch_errors,
# classifies each position as het good, het bad or non informative.
def determine_err_locs(positions, evaluated_positions, switches_indexes):

    position_switches = np.take(evaluated_positions, switches_indexes)
    result = np.zeros(positions.shape)
    #print position_switches
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
    hets = (calldata > 0).sum(axis=1)
    errs = (calldata == 2).sum(axis=1)
    bins = np.linspace(start, stop, int(pos.size/1000))

    midpoints = np.array([np.mean([x, y]) for x, y in zip(bins[:-1], bins[1:])])

    esites, _, _ = binned_statistic(pos, hets, np.sum, bins=bins)
    errors, _, _ = binned_statistic(pos, errs, np.sum, bins=bins)
    err_rate = errors/esites

    err_rate[esites < 10] = np.NaN

    fig = plt.figure(figsize=(16, 10))

    # plot genotypes
    ax = plt.subplot2grid((8, 1), (2, 0), rowspan=6)
    gt.plot_discrete_calldata(calldata, labels=[s.decode() for s in parent_ids],
                              colors=('grey', 'white', 'red'),
                              states=(0, 1, 2), ax=ax)

    ax = plt.subplot2grid((8, 1), (1, 0), rowspan=1)
    loc.plot_variant_locator(pos, step=pos.size/100, ax=ax, flip=True)
    ax.xaxis.set_ticklabels([])

    ax = plt.subplot2grid((8, 1), (0, 0), rowspan=1)
    ax.plot(midpoints, err_rate, 'r')
    ax.set_xlim((start, stop))
    ul = ceil(10 * np.nanmax(err_rate))/10
    ax.set_ylim((0, ul))
    ax.set_yticks(np.arange(0, ul + 0.01, 0.1))
    ax.set_ylabel("Switch error rate")

    ax.xaxis.set_major_formatter(
        ticker.FuncFormatter(lambda y, ps: '%.1f' % (y*1e-6)))
    fig.subplots_adjust(hspace=0.5)
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

parser.add_argument('--accessibility', '-A', action='store', default=None,
                    dest='accessibility', help='Accessibility h5', type=str,
                    required=False)

parser.add_argument('--samples', '-S', action='store', nargs="+", default=None,
                    dest='samples', required=False,
                    help='Which samples to evaluate.')

parser.add_argument('--chr', '-C', default=None, required=True,
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
    test_idx = [test_samples.index(s) for s in args.samples]
    eval_idx = [eval_samples.index(s) for s in args.samples]
    sample_names = args.samples
else:
    intersection = np.intersect1d(eval_samples, test_samples)
    test_idx = [test_samples.index(s) for s in intersection]
    eval_idx = [eval_samples.index(s) for s in intersection]
    sample_names = intersection

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

# equally keep only the eval positions that are in the test set
eval_keep_pos = np.in1d(e_pos, t_pos)
e_gt = e_gt.compress(eval_keep_pos, axis=0)
e_pos = np.compress(eval_keep_pos, e_pos)

assert e_gt.shape == t_gt.shape, ("Not same shape:", e_gt.shape, t_gt.shape)

is_missing = e_gt.is_missing()

data = dict()

for idx, sid in enumerate(sample_names):

    hz = gt.is_het(e_gt[:, idx])
    sample_pos = np.compress(hz, e_pos)
    sample_gt = np.compress(hz, t_gt[:, idx], axis=0)
    sample_gs = np.compress(hz, e_gt[:, idx], axis=0)

    switch = ph.switch.determine_haplotype_switch(sample_gs, sample_gt)
    pos_sw = ph.switch.derive_position_switch_array(switch)
    data[sid] = switch, pos_sw, sample_pos

het_positions = gt.is_het(e_gt).any(axis=1)
pl = plot_switch_errors(np.compress(het_positions, e_pos), data)
pl.savefig(filestem + "switch_errors.png")

# now build the windowed MSD array.
# ignore ROH for now....
#windows = np.arange(1, contig_lengths[args.chrom], args.winsize).astype('int')

# mmd, msd, conf, nmarker, nerror = calculate_switch_distances(
#     windows, data, window_size=5e5, roh=None)