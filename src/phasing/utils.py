# Collection of functions for phasing analysis/interpretation.

# external
import numpy as np
import matplotlib.pyplot as plt

# internal
import anhima
import sh
from scipy.special import gammaln
from itertools import izip
import pandas as pd
import hashlib


def get_relevant_haplotpes(pat_idx, mat_idx, pro_idx, hap_gt, hap_pos):

    """ for a given indices returns the parental and progeny genotypes
        cutting out any non segregating sites for plotting

     - pat_idx index of paternal diploid genotype
     - mat_idx index of maternal diploid genotype
     - pro_idx indices of progeny diploid genotypes
         (all indices refer to the 2nd dimension of hap_gt)
     - hap_gt 3D genotype matrix, variants X samples X ploidy.
     - hap_pos a numpy array of the positions, this is required so we can cut
       out the non interesting positions.

    - returns the maternal, paternal, progeny and position genotypes

    """
    # obtain haplotypes
    pat_haplotypes = np.compress(pat_idx, hap_gt, axis=1).reshape(
        hap_gt.shape[0], -1)
    mat_haplotypes = np.compress(mat_idx, hap_gt, axis=1).reshape(
        hap_gt.shape[0], -1)
    
    progeny_haps = np.compress(pro_idx, hap_gt, axis=1).reshape(
        hap_gt.shape[0], -1)

    # before returning, include only segregating positions
    is_seg = anhima.gt.is_het(pat_haplotypes) | anhima.gt.is_het(mat_haplotypes)
    
    pat_haplotypes = np.compress(is_seg, pat_haplotypes, axis=0)
    mat_haplotypes = np.compress(is_seg, mat_haplotypes, axis=0)
    progeny_haps = np.compress(is_seg, progeny_haps, axis=0)
    cmp_pos = np.compress(is_seg, hap_pos, axis=0)
        
    return pat_haplotypes, mat_haplotypes, progeny_haps, cmp_pos


def calc_distance_matrix(parental_haplotypes, progeny_haplotypes):

    """
:param parental_haplotypes: matrix of parental haplotypes
:param progeny_haplotypes: matrix of progeny haplotypes
:return: distance matrix from parental haplotypes. Just counts the number
 of non-agreements. Resulting matrix is S x 4 (4 parental haplotypes)
"""

    # we exclude homozygotic sites as they are uninformative.
    is_het = anhima.gt.is_het(parental_haplotypes)
    par_hts = np.compress(is_het, parental_haplotypes, axis=0)
    pro_hts = np.compress(is_het, progeny_haplotypes,  axis=0)

    # I suspect there is a far easier way to do this!
    distance_mat = [np.apply_along_axis(
        lambda a, b: np.absolute(a-b).sum(), 
        0, 
        pro_hts,
        par_hts[:, i]
    ) for i in range(par_hts.shape[1])]

    distance_mat = np.array(distance_mat).T

    # Normalize by the number of hets
    return distance_mat / float(is_het.sum())


def calc_best_assignment_qual(parental_haplotypes, progeny_haplotypes):

    """
function that calculates number of impossible haplotype assignments.
:param parental_haplotypes: 2D array, N x 4
:param progeny_haplotypes: 2D array of same N x S. Every pair is an individual.
:return: Returns the vector of how parimonious the best explantion is.
         Result is length S/2
"""

    distance_mat = calc_distance_matrix(parental_haplotypes, progeny_haplotypes)
    
    dmr = distance_mat.reshape(distance_mat.shape[0]/2, 2, 4)
    
    combs = [(0, 2), (0, 3), (1, 2), (1, 3)]
    
    l = []
    for i in range(dmr.shape[0]):
        v = dmr[i]
        l.append(np.min(
            [v[0, c[0]] + v[1, c[1]] for c in combs]
        ))
    assert len(l) == progeny_haplotypes.shape[1]/2

    return np.mean(l)


def plot_haplotype_descent(parental_haplotypes, progeny_haplotypes, positions):

    """
Creates a plot of haplotype descent using matplot lib
:param parental_haplotypes: 2D array. N x 2
:param progeny_haplotypes:  2D array. N x S
:param positions: 1D array of same N
:return: Nothing
"""

    n_parents = parental_haplotypes.shape[1]
    n_progeny = progeny_haplotypes.shape[1]
    
    # determine figure dimensions
    # plot haplotypes, allow 1 for each progeny, 1.5 for each parent
    # and 2 for pos
    fig_h = 1.0 + 0.3*(n_progeny + n_parents)
    par_pro_gap = 0.02
    par_par_gap = 0.01
    hap_width = (1.0-par_pro_gap)/(n_parents+n_progeny)

    pat_panel = (0, 0.1 + (2+n_progeny)*hap_width+par_par_gap+par_pro_gap,
                 1, 2*hap_width)
    mat_panel = (0, 0.1 + (n_progeny*hap_width)+par_pro_gap,
                 1, 2*hap_width)
    mid_panel = (0, 0.1, 1, n_progeny*hap_width)
    bot_panel = (0, 0, 1, 0.08)
    
    fig = plt.figure(figsize=(14, fig_h))

    ax = fig.add_axes(pat_panel)
    anhima.gt.plot_diploid_genotypes(
        np.flipud(parental_haplotypes[:, 0:2]),
        ax=ax, 
        colors=('white', 'black', 'pink')
    )
    
    ax = fig.add_axes(mat_panel)
    anhima.gt.plot_diploid_genotypes(
        np.flipud(parental_haplotypes[:, 2:4]),
        ax=ax, 
        colors=('white', 'black', 'pink')
    )
    
    ax = fig.add_axes(mid_panel)
    anhima.gt.plot_diploid_genotypes(
        np.flipud(progeny_haplotypes), 
        ax=ax, 
        colors=('white', 'black', 'pink')
    )

    # plot variant locator
    ax = fig.add_axes(bot_panel)
    anhima.loc.plot_variant_locator(positions, step=10, ax=ax)


def create_sh_script(filename, commands=None, outfile=None):

    """
# Not strictly a phasing function, but so commonly used, may as well put here!
:param filename: name of file to write
:param commands: list of strings that will be system executed
:param outfile: optional, creates an 'ok' version of this file
:return: None
"""

    # create script file
    if outfile is None:
        outfile = filename
    
    touch_cmd = "touch {FILE}"
    
    script = open(filename, 'w')
    script.write("#! /bin/bash" + "\n")
    script.write("set -e" + "\n")
    script.write("set -o pipefail" + "\n")
    for cmd in commands:
        script.write(cmd + "\n")
    script.write(touch_cmd.format(FILE=outfile + '.ok') + "\n")
    sh.chmod('+x', filename)
    script.close()


def calc_regions(size, nbins=20, overlap=0):

    """
    :param size: size of total region
    :param nbins: number of regions
    :param overlap: stagger
    :return:
    """
    if overlap is None:
        overlap = size/(nbins*20)
    approx_size = 1 + size/nbins + overlap - (overlap/nbins)
    print approx_size, overlap
    regions = []

    for i in range(nbins):
        start = i*(approx_size-overlap) + 1
        stop = start + approx_size
        regions.append((start, stop))
    return regions


def determine_switches(a):
    return np.diff(np.where(np.concatenate(([a[0]], a[:-1] != a[1:],
                                            [True])))[0])


# used to remove single errors
def forgive(a, ignore=1):
    assert 2 > ignore >= 0
    idx = np.where(np.concatenate(([a[0]], a[:-1] != a[1:], [True])))[0]
    switch = np.diff(idx)
    forgive_me = np.ones(a.shape, dtype='bool')
    forgive_me[idx[switch <= ignore]] = False
    return forgive_me


def calculate_switch_error(inheritance, ignore_size=0):

    # only 1s and 2s are relevant
    exclude = np.any(inheritance < 3, axis=1)
    inh_copy = np.compress(exclude, inheritance.copy(), axis=0)

    forgiven = [forgive(col, ignore_size) for col in inh_copy.T]
    switches = [determine_switches(np.compress(fgv, col))
                for col, fgv in zip(inh_copy.T, forgiven)]

    switch_e = [s.size - 1 for s in switches]
    ignored = [np.sum(~f) for f in forgiven]

    return switch_e, ignored, inh_copy.shape


def calculate_switch_length(inheritance, positions, ignore_size=0,
                            index_only=False):

    # only 1s and 2s are relevant
    exclude = np.any(inheritance < 3, axis=1)
    inh_copy = np.compress(exclude, inheritance.copy(), axis=0)

    forgiven = [forgive(col, ignore_size) for col in inh_copy.T]
    switches = [determine_switches(np.compress(fgv, col))
                for col, fgv in zip(inh_copy.T, forgiven)]

    if index_only:
        mean_length = [np.mean(s) for s in switches]
        medi_length = [np.median(s) for s in switches]
        maxi_length = [np.median(s) for s in switches]
    else:
        assert inheritance.shape[0] == positions.shape[0]
        pos = np.compress(exclude, positions)

        filtered_pos = [np.insert(np.take(np.compress(fgv, pos),
                                          sw.cumsum() - 1), 0, pos[0])
                        for fgv, sw in zip(forgiven, switches)]

        mean_length = np.array([np.mean(np.diff(f)) for f in filtered_pos])
        medi_length = np.array([np.median(np.diff(f)) for f in filtered_pos])
        maxi_length = np.array([np.max(np.diff(f)) for f in filtered_pos])

    return mean_length, medi_length, maxi_length


def plot_ped_haplotype_inheritance(parent_genotypes,
                                   progeny_genotypes,
                                   positions,
                                   filename=None,
                                   downsample=10000,
                                   select='eitherhet',
                                   title='Haplotype inheritance',
                                   inheritance_colors=('red', 'blue',
                                                       'green', 'orange',
                                                       'black', 'yellow',
                                                       'white'),
                                   spacer=0.05,
                                   panel_height_ratios=(3.0, 3.0, 1.0, 1.0),
                                   progeny_labels=None):
    """
    Creates a plot for each pedigree in the pedigree dict
    :param run: an instance of Tool
    :param pedigree: a dict
    :param inheritance_colors: colour scheme
    :param spacer: gap between panels
    :param panel_height_ratios: relative size of each panel
    :return:
    """

    panel_heights = np.array(panel_height_ratios)/np.sum(panel_height_ratios)\
        - spacer

    selection = np.ones(parent_genotypes.shape[0], dtype='bool')
    if select == 'eitherhet':
        selection = anhima.gt.is_het(parent_genotypes).any(axis=1)
    elif select == 'bothhet':
        selection = anhima.gt.is_het(parent_genotypes).all(axis=1)

    # downsample
    if downsample is not None:
        where = np.where(selection)[0]
        idx = np.random.choice(where, downsample, replace=False)
        toplot = np.zeros(selection.shape, dtype='bool')
        toplot[idx] = True
        print 'Downsampling from {0} to {1} for plotting'.format(str(
            selection.sum()), str(downsample))
    else:
        toplot = selection

    father = np.compress(toplot, parent_genotypes[:, 0], axis=0)
    mother = np.compress(toplot, parent_genotypes[:, 1], axis=0)

    # NB: critical assumption of genotypes here
    paternal_haplotypes = np.compress(toplot, progeny_genotypes,
                                      axis=0)[:, :, 0]
    maternal_haplotypes = np.compress(toplot, progeny_genotypes,
                                      axis=0)[:, :, 1]

    positions = np.compress(selection, positions, axis=0)

    # Again assumption of correct location of parents in dict
    maternal_inheritance = anhima.ped.diploid_inheritance(
        mother, maternal_haplotypes)
    paternal_inheritance = anhima.ped.diploid_inheritance(
        father, paternal_haplotypes)

    axes = [(0, n*spacer + panel_heights[:n].sum(),
             1, panel_heights[n]) for n in range(panel_heights.size)]

    fig, ax = plt.subplots(figsize=(12, 8))
    progeny_labels.reverse()

    plt.title(title)   # subplot 211 title

    # (left, bottom, width, height)
    ax = fig.add_axes(axes.pop())
    anhima.loc.plot_windowed_variant_density(positions,
                                             window_size=50000,
                                             ax=ax)

    ax = fig.add_axes(axes.pop())
    anhima.loc.plot_variant_locator(positions,
                                    step=1000,
                                    ax=ax,
                                    flip=False)

    ax = fig.add_axes(axes.pop())
    anhima.gt.plot_discrete_calldata(paternal_inheritance,
                                     colors=inheritance_colors,
                                     labels=progeny_labels,
                                     states=range(1, 8),
                                     ax=ax)

    ax = fig.add_axes(axes.pop())
    anhima.gt.plot_discrete_calldata(maternal_inheritance,
                                     colors=inheritance_colors,
                                     labels=progeny_labels,
                                     states=range(1, 8),
                                     ax=ax)
    if filename is not None:
        plt.savefig(filename, bbox_inches='tight')
    return ax


def plot_single_hap_inheritance(parent_genotypes, gamete_haplotypes, positions,
                                show_all=False, filename=None, downsample=10000,
                                inheritance_colors=('red', 'blue', 'green',
                                                    'orange', 'black', 'yellow',
                                                    'white'),
                                spacer=0.05,
                                phr=(1.0, 1.0, 4.0),
                                progeny_labels=None):
    """
    Creates a plot for each pedigree in the pedigree dict
    :param parent_genotypes: an n x 2 array or a n x 1x 2 array that is squeez
    :param gamete_haplotypes: an n x S x 2 array
    :param inheritance_colors: colour scheme
    :param spacer: gap between panels
    :param phr: relative size of each panel
    :return:
    """

    panel_heights = np.array(phr)/np.sum(phr) - spacer

    selection = np.ones(parent_genotypes.shape[0], dtype='bool')
    if not show_all:
        selection = anhima.gt.is_het(parent_genotypes)

    # downsample
    if downsample is not None:
        where = np.where(selection)[0]
        idx = np.random.choice(where, downsample, replace=False)
        toplot = np.zeros(selection.shape, dtype='bool')
        toplot[idx] = True
        print 'Downsampling from {0} to {1} for plotting'.format(str(
            selection.sum()), str(downsample))
    else:
        toplot = selection

    print "Plotting {0} variants".format(toplot.sum())

    parent = np.compress(toplot, parent_genotypes, axis=0)

    # NB: critical assumption of genotypes here
    haplotypes = np.compress(toplot, gamete_haplotypes, axis=0)
    positions = np.compress(selection, positions, axis=0)

    # Again assumption of correct location of parents in dict
    inheritance = anhima.ped.diploid_inheritance(parent, haplotypes)

    axes = [(0, n*spacer + np.sum(panel_heights[:n]),
             1, panel_heights[n]) for n in range(np.size(panel_heights))]

    fig, ax = plt.subplots(figsize=(12, 8))

    if progeny_labels is None:
        progeny_labels = [''] * gamete_haplotypes.shape[1]
    else:
        progeny_labels.reverse()

    ax = fig.add_axes(axes.pop())
    anhima.gt.plot_discrete_calldata(inheritance,
                                     colors=inheritance_colors,
                                     labels=progeny_labels,
                                     states=range(1, 8),
                                     ax=ax,
                                     pcolormesh_kwargs={'edgecolors': 'black'})

    ax = fig.add_axes(axes.pop())
    try:
        anhima.loc.plot_variant_locator(positions,
                                        step=toplot.sum()/100,
                                        ax=ax,
                                        flip=False)
    except ValueError:
        pass

    # (left, bottom, width, height)
    ax = fig.add_axes(axes.pop())
    try:
        window = (positions[-1] - positions[0])/100
        anhima.loc.plot_windowed_variant_density(positions,
                                                 window_size=window,
                                                 ax=ax)
    except ValueError:
        pass

    if filename is not None:
        plt.savefig(filename, bbox_inches='tight')
    return ax


def read_pedigree_table(path, pedigree_id_col='cross', status_id_col='role',
                        sep="\t"):
    # This file defines the known pedigrees we are about to test.
    # expecting data in the format of:
    # id<TAB>cross<TAB>role
    # A1<TAB>19-2<TAB>parent
    ped_tbl = pd.read_csv(path, sep=sep, index_col=1)

    pedigree = dict()
    for x in set(ped_tbl[pedigree_id_col]):
        pedigree[x] = {k: [] for k in set(ped_tbl[status_id_col])}

    # loop through each row of table
    for sid, row in ped_tbl.iterrows():
        pedigree[row[pedigree_id_col]][row[status_id_col]].append(sid)

    return pedigree, ped_tbl


def md5_for_file(f, block_size=2**20):
    fh = open(f, 'rb')
    md5 = hashlib.md5()
    while True:
        data = fh.read(block_size)
        if not data:
            break
        md5.update(data)
    return md5.hexdigest()


def create_samples_file(path=None, output=None, samples=None):
    pedigree, ped_tbl = read_pedigree_table(path)
    fh = open(output, mode='w')
    header_1 = " ".join(['ID_1', 'ID_2', 'missing', 'father',
                         'mother', 'sex', 'plink_pheno'])
    header_2 = " ".join(['0', '0', '0', 'D', 'D', 'D', 'B'])

    if samples is None:
        samples = ped_tbl.index.tolist()

    fh.write(header_1 + "\n")
    fh.write(header_2 + "\n")

    ped_counter = {p: 0 for p in set(ped_tbl.cross.tolist() + ['none'])}

    for s in samples:

        maternal_id = '0'
        paternal_id = '0'
        cross = 'none'

        try:
            match = ped_tbl.loc[s]
            cross = match['cross']
            parents = ped_tbl[(ped_tbl.cross == cross) & (ped_tbl.role ==
                                                          'parent')].index
            if match.name not in parents:
                paternal_id = parents[0]
                maternal_id = parents[1]

        except KeyError:
            pass

        ped_counter[cross] += 1
        family = cross + '_' + str(ped_counter[cross])

        line = " ".join([family, s, '0', paternal_id, maternal_id, '0', '-9'])
        fh.write(line + "\n")


def return_classification(geno_data):
    return np.prod(geno_data + 1, axis=1)


def log_factorial(x):
    """Returns the logarithm of x!
    Also accepts lists and NumPy arrays in place of x."""
    return gammaln(np.array(x)+1)


def memoize(function):
    memo = {}

    def wrapper(*inargs):
        args = np.vstack(inargs).tostring()
        if args in memo:
            return memo[args]
        else:
            rv = function(*inargs)
            memo[args] = rv
            return rv
    return wrapper


@memoize
def log_multinomial(xs, ps):
    xs, ps = np.array(xs), np.array(ps)
    assert ps.sum() == 1
    result = log_factorial(np.sum(xs)) - np.sum(log_factorial(xs)) + \
        np.sum(xs * np.log(ps))
    return result


def get_error_likelihood(parental_genotypes, progeny_genotypes, pe=0.001):
    """Returns the log likelihood that the stated parental genotypes are
     correct. """

    # classification = {1: 'HomRef_HomRef', 2: 'HomRef_Het',
    #                   3: 'HomRef_HomAlt', 4: 'Het_Het',
    #                   6: 'HomAlt_Het', 9: 'HomAlt_HomAlt'}

    lookup = {1: (1-2*pe, pe, pe),
              2: (0.5-pe/2, 0.5-pe/2, pe),
              3: (pe, 1-2*pe, pe),
              4: (0.25, 0.5, 0.25),
              6: (pe, 0.5-pe/2, 0.5-pe/2),
              9: (pe, pe, 1-2*pe)}

    parental_genotypes = anhima.gt.as_012(parental_genotypes)
    progeny_genotypes = anhima.gt.as_012(progeny_genotypes)

    counts = np.vstack([np.sum(progeny_genotypes == 0, axis=1),
                        np.sum(progeny_genotypes == 1, axis=1),
                        np.sum(progeny_genotypes == 2, axis=1)])

    classification = return_classification(parental_genotypes)
    assert classification.ndim == 1
    res = list()
    for i in xrange(classification.size):
        if classification[i] == 0:
            res.append(0.0)
            continue
        r = log_multinomial(counts[:, i],
                            lookup[classification[i]])

        v = np.max([log_multinomial(counts[:, i], lookup[key]) for key
                    in lookup.keys() if key != classification[i]])
        res.append(r - v)
    return np.array(res)


def get_consecutive_true(a):
    if a.sum() == 0:
        return 0
    else:
        return np.diff(np.where(np.concatenate(([a[0]],
                                                a[:-1] != a[1:],
                                                [True])))[0])[::2].max()


def mask_2d(a, c, threshold=0, value=False):
    mask = np.array(c < threshold)
    cp_a = a.copy()
    cp_a[mask] = value
    return cp_a, mask


def mask_3d(a, c, threshold=0, value=(-1, -1)):
    mask = np.array(c < threshold)
    cp_a = a.copy()
    cp_a[mask] = value
    return cp_a, mask