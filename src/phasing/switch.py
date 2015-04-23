__author__ = 'Nicholas Harding'

import numpy as np
from anhima.gt import is_het, is_hom_alt, is_hom_ref
#from anhima.ped import diploid_inheritance
#from anhima.ped import diploid_mendelian_error_biallelic

# collection of functions that allow calculation of switch error rates etc
# this is under active development!!


def check_equivalency(a, b):
    assert a.ndim == b.ndim == 3
    assert np.array_equal(a.shape, b.shape)
    assert np.array_equal(a.sum(axis=2), b.sum(axis=2))


def determine_haplotype_switch(geno_a, geno_b):
    """
    Operates on two alternative haplotype assertions on two identical genotypes.
    Returns haplotype switch array, as integers to avoid confusion
    The information is encoded in runs of 0s and 1s. a 0 indicates that the
    0th allele from the first genotype array matches the 0th allele from the
    second. A 1 indicates that the 0th allele from the first matches the
    first from the second. As the order of the 3rd axis is arbitrary the
    distinction between 0/1s is unimportant, we may as well work from the 1st
    axis and invert the result. Will error if any inconsistencies.

    :return: np.array() of length N
    """

    # check for equivalence
    check_equivalency(geno_a, geno_b)
    return compute_allele_agreement(np.take(geno_a, 0, axis=2), geno_b)


def compute_allele_agreement(hap, geno):

    a = hap == np.take(geno, 0, axis=2)
    b = hap == np.take(geno, 1, axis=2)
    # a = 0, b = 1, a & b == err, not a | b = 2.

    out = 4 - ((a + 1) * (b + 2))
    # 2 = unkn, 1 = 2nd allele, 0 = 1st allele, -2 = Fail- both!
    assert not np.any(-2 == out), "Error: some homozygotes must be present"
    return out


def determine_inherited_allele(hap, geno, allow_mendelian_errors=False):
    """
    Similar to above, except takes a genotype and a haplotype instead of 2
    genotypes. Not putatively identical genotypes. If the allow flag is True,
    then errors are NOT marked as switches. but reported as errors separately
    :param geno:
    :param hap:
    :return:
    """

    # check there are no inconsistencies ie 1 from 0/0.
    mendel_incon = (is_hom_ref(geno) & (1 == hap)) | \
                   (is_hom_alt(geno) & (0 == hap))

    if not allow_mendelian_errors:
        assert not mendel_incon.any(), "Mendelian inconsistencies observed"

    agreement = compute_allele_agreement(hap, geno)

    assert np.array_equal(agreement == 2, mendel_incon)
    return agreement


def derive_position_switch_array(alleles):
    # The position switch array describes the runs of allele ids.
    # Any non inherited alleles are ignored...just masked.
    # consistent_alleles = alleles != 2
    # a = np.compress(consistent_alleles, alleles)
    # if consistent_alleles.any():
    #     print "{0} non consistent site(s), either mutation or genotype " \
    #           "error. These are ignored for purposes of switches.".format(
    #           np.invert(consistent_alleles).sum())
    #
    m = alleles[:-1] != alleles[1:]
    co = np.where(np.concatenate(([False], m, [True])))[0]
    ans = np.diff(np.insert(co, 0, 0))
    return ans


def adjust_position_switch_array(pos_sw, tolerance=1):

    a = pos_sw.copy()
    while True:

        loci = np.diff(pos_sw) <= tolerance
        print "loci:", loci
        if not loci.any():
            break
        xx = np.where(loci)[0][0]
        print "xx:", xx
        a[xx] = a[xx - 1]
    return pos_sw


def count_switches(pos_sw):
    return pos_sw.size - 1


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
    switches = [derive_position_switch_array(np.compress(fgv, col))
                for col, fgv in zip(inh_copy.T, forgiven)]

    switch_e = np.array([s.size - 1 for s in switches])
    ignored = np.array([np.sum(~f) for f in forgiven])

    return switch_e, ignored, inh_copy.shape, switches


def calculate_switch_length(inheritance, positions, ignore_size=0,
                            index_only=False):
    assert inheritance.shape[0] == positions.size

    # only 1s and 2s are relevant
    exclude = np.any(inheritance < 3, axis=1)
    inh_copy = np.compress(exclude, inheritance.copy(), axis=0)

    forgiven = [forgive(col, ignore_size) for col in inh_copy.T]
    switches = [derive_position_switch_array(np.compress(fgv, col))
                for col, fgv in zip(inh_copy.T, forgiven)]

    filtered_pos = None
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

    return mean_length, medi_length, maxi_length, filtered_pos

