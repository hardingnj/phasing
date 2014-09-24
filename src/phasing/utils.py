# Collection of functions for phasing analysis/interpretation.

# external
import numpy as np
import matplotlib.pyplot as plt

# internal
import anhima
import sh


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


def calc_windows(k, end, contig, start=1):
    """
    # Define function that detemines the regions, in chunks to split the jobs
    # Function returns a list of strings
    :param k: desired size
    :param end: final base (inclusive)
    :param contig: the contig
    :param start: start of the window, usually 0/1
    :return:
    """

    regions = []
    i = start
    while i < end:
        regions.append(contig + '{S1}'
                       + str(i) + "{S2}"
                       + str(min(end, i + k - 1)))
        i += k

    return regions