# Collection of functions for phasing analysis/interpretation.

# external
import numpy as np
import matplotlib.pyplot as plt

# internal
import anhima.gt
import anhima.loc

# function that for a given cross returns the par and prog genotypes
# - pat_idx index of paternal diploid genotype
# - mat_idx index of maternal diploid genotype
# - pro_idx indices of progeny diploid genotypes
#     (all indices refer to the 2nd dimension of hap_gt)
# - hap_gt 3D genotype matrix, variants X samples X ploidy.
# - hap_pos a numpy array of the positions, this is required so we can cut out the non interesting positions.

def get_relevant_haplotpes(pat_idx, mat_idx, pro_idx, hap_gt, hap_pos):
    
    # obtain haplotypes
    pat_haplotypes  = np.compress(pat_idx, hap_gt, axis=1).reshape(hap_gt.shape[0],-1)
    mat_haplotypes  = np.compress(mat_idx, hap_gt, axis=1).reshape(hap_gt.shape[0],-1)
    
    progeny_haps = np.compress(pro_idx, hap_gt, axis=1).reshape(hap_gt.shape[0],-1)

    # before returning, include only segregating positions
    is_seg = anhima.gt.is_het(pat_haplotypes) | anhima.gt.is_het(mat_haplotypes)
    
    pat_haplotypes = np.compress(is_seg, pat_haplotypes, axis = 0)
    mat_haplotypes = np.compress(is_seg, mat_haplotypes, axis = 0)
    progeny_haps   = np.compress(is_seg, progeny_haps, axis = 0)
    cmp_pos        = np.compress(is_seg, hap_pos, axis = 0)
        
    return pat_haplotypes, mat_haplotypes, progeny_haps, cmp_pos

# -----------------

# CALC_MEAN_EUCLIDEAN_DIST
# function that calculates a distance matrix from parental haplotypes. 
#  Just counts the number of non-agreements
# - parental_haplotypes a 2D array, N x 4
# - progeny_haplotypes a 2D array of same N x S
# Resulting matrix is S x 4 (4 parental haplotypes)

def calc_distance_matrix(parental_haplotypes, progeny_haplotypes):

    # we exclude homozygotic sites as they are uninformative.
    is_het = anhima.gt.is_het(parental_haplotypes)
    par_hts = np.compress(is_het, parental_haplotypes, axis = 0)
    pro_hts = np.compress(is_het, progeny_haplotypes,  axis = 0) 

    distance_mat = [np.apply_along_axis(
        lambda a, b: np.absolute(a-b).sum(), 
        0, 
        pro_hts,
        par_hts[:,i]
    ) for i in range(par_hts.shape[1])]

    distance_mat = np.array(distance_mat).T
    return distance_mat / float(is_het.sum())

# -----------------

# CALC_BEST_ASSIGNEMENT_QUAL
# function that calculates number of impossible haplotype assignments.
# - parental_haplotypes a 2D array, N x 4
# - progeny_haplotypes a 2D array of same N x S. Every pair is an individual.

# Returns the vector of how parimonious the best explantion is. Result is length S/2 

def calc_best_assignment_qual(parental_haplotypes, progeny_haplotypes):

    distance_mat = calc_distance_matrix(parental_haplotypes, progeny_haplotypes) 
    
    dmr = distance_mat.reshape(distance_mat.shape[0]/2, 2, 4)
    
    combs = [(0,2), (0,3), (1,2), (1,3)]
    
    l = []
    for i in range(dmr.shape[0]):
        v = dmr[i]
        l.append(np.min(
            [v[0, c[0]] + v[1, c[1]] for c in combs]
        ))
    assert len(l) == progeny_haplotypes.shape[1]/2

    return np.mean(l)

# ------------

# Plot haplotypes
# Creates a plot of haplotype descent using matplot lib

def plot_haplotype_descent(parental_haplotypes, progeny_haplotypes, positions):

    n_parents = parental_haplotypes.shape[1]
    n_progeny = progeny_haplotypes.shape[1]
    
    # determine figure dimensions
    # plot haplotypes, allow 1 for each progeny, 1.5 for each parent and 2 for pos
    fig_h = 1.0 + 0.3*(n_progeny + n_parents)
    par_pro_gap = 0.02
    par_par_gap = 0.01
    hap_width = (1.0-par_pro_gap)/(n_parents+n_progeny)

    pat_panel = (0, 0.1 + (2+n_progeny)*hap_width+par_par_gap+par_pro_gap, 1, 2*hap_width)
    mat_panel = (0, 0.1 + (n_progeny*hap_width)+par_pro_gap,   1, 2*hap_width)
    mid_panel = (0, 0.1,                                       1, n_progeny*hap_width)
    bot_panel = (0, 0,                                         1, 0.08)
    
    fig = plt.figure(figsize=(14, fig_h))

    ax = fig.add_axes(pat_panel)
    anhima.gt.plot_diploid_genotypes(
        np.flipud(parental_haplotypes[:,0:2]), 
        ax=ax, 
        colors = ('white', 'black', 'pink')
    )
    
    ax = fig.add_axes(mat_panel)
    anhima.gt.plot_diploid_genotypes(
        np.flipud(parental_haplotypes[:,2:4]), 
        ax=ax, 
        colors = ('white', 'black', 'pink')
    )
    
    ax = fig.add_axes(mid_panel)
    anhima.gt.plot_diploid_genotypes(
        np.flipud(progeny_haplotypes), 
        ax=ax, 
        colors = ('white', 'black', 'pink')
    )

    # plot variant locator
    ax = fig.add_axes(bot_panel)
    anhima.loc.plot_variant_locator(positions, step=10, ax=ax);


