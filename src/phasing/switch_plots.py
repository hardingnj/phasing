__author__ = 'Nicholas Harding'

import numpy as np
import anhima
import matplotlib.pyplot as plt


def plot_hap_switches(block_tbl, het_positions, index_positions, filename=None,
                      title=None, color=('grey', 'red', 'blue')):

    """
    Creates a plot for each pedigree in the pedigree dict
    :param block_tbl: output of anhima's tabulate_block
    :param het_positions: n x 1 array with each het position in the parent
    :param color: colour scheme
    :return: axis of created plot
    """

    prog = block_tbl.groupby(level=0)
    n_samples = len(prog.groups.keys())
    contig_len = block_tbl.stop_max_pos.max()
    index_limit = block_tbl.stop_max_idx.max()

    fig = plt.figure(figsize=(14, int(n_samples * 0.3) + 5))

    for i, p in enumerate(prog.groups.keys()):
        ax = plt.subplot2grid((n_samples + 5, 1), (i, 0), rowspan=1)
        ax.xaxis.set_ticklabels([])
        ax.set_xlim(0, index_limit)
        ax.set_ylim(0, 1)
        ax.set_yticks([0.5])
        ax.yaxis.set_tick_params(size=0)
        ax.yaxis.set_ticklabels([p])

        x = prog.get_group(p)
        for j, row in x.iterrows():
            low = row.start_max_idx
            upp = row.stop_min_idx
            _ = ax.axvspan(low, upp, facecolor=color[row.state],
                           alpha=0.7, edgecolor='none')

    ax.get_xaxis().set_visible(True)

    ax = plt.subplot2grid((n_samples + 5, 1), (n_samples + 1, 0), rowspan=2)
    ax.set_xticks([])
    anhima.loc.plot_variant_locator(index_positions, 1000, ax=ax, flip=False)

    ax = plt.subplot2grid((n_samples + 5, 1), (n_samples + 3, 0), rowspan=2)
    anhima.loc.plot_windowed_variant_density(het_positions, 1e5, ax=ax)
    _ = ax.set_xticks(np.arange(0, contig_len, 5e6))
    _ = ax.set_xticklabels(["{0:.1f}Mb".format(x/1e6)
                            for x in np.arange(0, contig_len, 5e6)])

    yticks = ax.yaxis.get_major_ticks()
    for i in range(1, len(yticks) - 1):
        yticks[i].label1.set_visible(False)

    if title is not None:
        fig.suptitle(title, fontsize=18)

    if filename is not None:
        plt.savefig(filename)

    return ax