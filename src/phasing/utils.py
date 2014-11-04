# Collection of functions for phasing analysis/interpretation.

# external
import numpy as np
import matplotlib.pyplot as plt

# internal
import anhima
import sh
import os
import re
import yaml
import algorithms
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


def reconstruct_run(directory):
    """
    This function returns the correct tool object with appropriate settings
    given an input directory
    :param directory:
    :return: an object of the correct type
    """

    # basename is run id
    directory = os.path.realpath(directory)
    assert os.path.isdir(directory)
    base, run_id = os.path.split(directory)

    param_yaml = os.path.join(directory, run_id + '_parameters.yaml')
    from_yaml = yaml.load(stream=open(param_yaml, 'r'))

    cls = getattr(algorithms, from_yaml['name'])

    parameters = from_yaml.get('parameters', [])

    # for backwards compatibility: grab values where
    if len(parameters) == 0:
        cl = re.compile('^-')
        from_yaml.pop('--output-max', None)
        for i, j in from_yaml.iteritems():
            if cl.match(i):
                parameters = parameters + [i, j]

    return cls(parameters=parameters,
               executable=from_yaml['executable'],
               outdir=from_yaml['base_dir'],
               version=from_yaml['version'],
               run_id=run_id)


# or genotypes/samples can be gathered from the run id. Also allows cacheing!
def calculate_pedigree_switch_error(run, pedigree, use_cache=True, **kwargs):
    """
    This function returns a tuple of panda dfs for mean and sd for each cross
    :param run: a Tool object, containing a parse_output method
    :param pedigree: a dict of pedigrees, with parent ids and progeny ids
    :param use_cache: whether to used cached results
    :param **kwargs passed to run.parse_output
    :return:
    """

    mean_fn = os.path.join('/tmp', run.run_id + 'switcherror_mean.csv')
    sd_fn = os.path.join('/tmp', run.run_id + 'switcherror_sd.csv')
    if os.path.exists(mean_fn) and os.path.exists(sd_fn) and use_cache:
        return pd.read_csv(mean_fn, index_col=0), \
            pd.read_csv(sd_fn, index_col=0)

    genotypes, samples, dic = run.parse_output(**kwargs)

    # output: 2 x pd dataframe
    # to do: make so in melted form: ie mat/pat as a indicator col.
    columns = ['cross', 'maternal', 'paternal']
    sd_df = pd.DataFrame(columns=columns)
    mean_df = pd.DataFrame(columns=columns)
    # row per cross, 2 columns, mat/pat, mean + sd
    for cross in pedigree.keys():
        parents = np.array([samples.index(p) for p in
                            pedigree[cross]['parent']])
        progeny = np.array([samples.index(p) for p in
                            pedigree[cross]['progeny']])

        # NB: critical assumption of genotypes here
        maternal_haplotypes = genotypes[:, progeny, 0]
        paternal_haplotypes = genotypes[:, progeny, 1]

        # Again assumption of correct location of parents in dict
        maternal_inheritance = anhima.ped.diploid_inheritance(
            genotypes[:, parents[0]], maternal_haplotypes)
        paternal_inheritance = anhima.ped.diploid_inheritance(
            genotypes[:, parents[1]], paternal_haplotypes)

        m_switcherror = calculate_switch_error(maternal_inheritance)
        p_switcherror = calculate_switch_error(paternal_inheritance)

        mean_df.loc[len(mean_df.index)] = (cross,
                                           m_switcherror.mean(),
                                           p_switcherror.mean())
        sd_df.loc[len(sd_df.index)] = (cross,
                                       m_switcherror.std(),
                                       p_switcherror.std())

    mean_df.to_csv(mean_fn)
    sd_df.to_csv(sd_fn)
    return mean_df, sd_df


def calculate_switch_error(inheritance):

    # only 1s and 2s are relevant
    exclude = np.any(inheritance < 3, axis=1)
    inheritance = np.compress(exclude, inheritance, axis=0)

    switch_sums = map(lambda (x, y): np.array(x != y, dtype='int8'),
                      izip(inheritance[1:], inheritance[0:-1]))

    return np.array(switch_sums).mean(axis=0)


def plot_ped_haplotype_inheritance(run, pedigree,
                                   inheritance_colors=('red', 'blue',
                                                       'green', 'orange',
                                                       'black', 'yellow',
                                                       'white'),
                                   spacer=0.02,
                                   panel_height_ratios=(3.0, 3.0, 1.0, 1.0)):
    """
    Creates a plot for each pedigree in the pedigree dict
    :param run: an instance of Tool
    :param pedigree: a dict
    :param inheritance_colors: colour scheme
    :param spacer: gap between panels
    :param panel_height_ratios: relative size of each panel
    :return:
    """

    genotypes, samples, dic = run.parse_output()
    panel_heights = np.array(panel_height_ratios)/np.sum(panel_height_ratios)\
        - spacer

    out_dir = os.path.join(run.outdir, 'plots')
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)

    for cross in pedigree.keys():
        # skip if fn exists
        filename = os.path.join(out_dir, cross + '.png')
        if os.path.isfile(filename):
            continue

        parents = np.array([samples.index(p) for p in
                            pedigree[cross]['parent']])
        progeny = np.array([samples.index(p) for p in
                            pedigree[cross]['progeny']])

        is_het = anhima.gt.is_het(genotypes[:, parents[0]]) | \
            anhima.gt.is_het(genotypes[:, parents[1]])

        mother = np.compress(is_het, genotypes[:, parents[0]], axis=0)
        father = np.compress(is_het, genotypes[:, parents[1]], axis=0)

        # NB: critical assumption of genotypes here
        maternal_haplotypes = np.compress(is_het,
                                          genotypes[:, progeny, 0],
                                          axis=0)
        paternal_haplotypes = np.compress(is_het,
                                          genotypes[:, progeny, 1],
                                          axis=0)

        positions = np.compress(is_het, dic['pos'], axis=0)

        # Again assumption of correct location of parents in dict
        maternal_inheritance = anhima.ped.diploid_inheritance(
            mother, maternal_haplotypes)
        paternal_inheritance = anhima.ped.diploid_inheritance(
            father, paternal_haplotypes)

        axes = [(0, n*spacer + panel_heights[:n].sum(),
                 1, panel_heights[n]) for n in range(panel_heights.size)]

        fig, ax = plt.subplots(figsize=(12, 8))

        progeny_labels = [p for p in pedigree[cross]['progeny'] if p in samples]
        progeny_labels.reverse()

        # (left, bottom, width, height)
        ax = fig.add_axes(axes.pop())
        anhima.loc.plot_windowed_variant_density(positions,
                                                 window_size=50000,
                                                 ax=ax)

        ax = fig.add_axes(axes.pop())
        anhima.loc.plot_variant_locator(positions,
                                        step=1000,
                                        ax=ax,
                                        flip=True)

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

        plt.savefig(filename, bbox_inches='tight')
        plt.close()


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