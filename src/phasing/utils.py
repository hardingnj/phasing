# Collection of functions for phasing analysis/interpretation.

# external
import numpy as np
import matplotlib.pyplot as plt
import gzip

# internal
import anhima
import sh
from scipy.special import gammaln
from collections import Counter
import pandas as pd
import hashlib
is_het = anhima.gt.is_het


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


# may get moved to py_shapeit
def parse_duohmm_genotype_error(fn, sample_names, threshold=0.5):

    samples = list()
    with gzip.open(fn, 'r') as fh:
        header = fh.readline().rstrip().split()

        for li in fh.readlines():
            ID, father, mother, RSID, pos, prob = li.rstrip().split()
            if float(prob) >= threshold:
                samples.append(ID)

    xx = Counter(samples)
    return np.array(map(xx.get, sample_names))


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


def log_multinomial(xs, ps):
    xs, ps = np.array(xs), np.array(ps)
    assert ps.sum() == 1.0
    result = log_factorial(np.sum(xs)) - np.sum(log_factorial(xs)) + \
        np.sum(xs * np.log(ps))
    return result


def get_error_likelihood(parental_genotypes, progeny_genotypes,
                         lookup=None, pe=0.001):
    """Returns the log likelihood that the stated parental genotypes are
     correct. """

    # classification = {1: 'HomRef_HomRef', 2: 'HomRef_Het',
    #                   3: 'HomRef_HomAlt', 4: 'Het_Het',
    #                   6: 'HomAlt_Het', 9: 'HomAlt_HomAlt'}
    if lookup is None:
        lookup = {1: (1-2*pe, pe, pe),
                  2: (0.5-pe/2, 0.5-pe/2, pe),
                  3: (pe, 1-2*pe, pe),
                  4: (0.25, 0.5, 0.25),
                  6: (pe, 0.5-pe/2, 0.5-pe/2),
                  9: (pe, pe, 1-2*pe)}

    parental_genotypes = anhima.gt.as_012(parental_genotypes)
    progeny_genotypes = anhima.gt.as_012(progeny_genotypes)

    memo_hash = dict()

    counts = np.vstack([np.sum(progeny_genotypes == 0, axis=1),
                        np.sum(progeny_genotypes == 1, axis=1),
                        np.sum(progeny_genotypes == 2, axis=1)]).T

    classification = return_classification(parental_genotypes)
    assert classification.ndim == 1
    res = np.zeros(classification.size)
    for i in xrange(classification.size):
        if classification[i] == 0:
            continue

        key = "_".join(np.insert(
            counts[i], 0, classification[i]).astype('string'))

        if not key in memo_hash:
            r = log_multinomial(counts[i], lookup[classification[i]])

            # now we take the keys from the other groups:
            alt_list = list()
            for pgt in filter(lambda k: k != classification[i], lookup.keys()):
                v_key = "_".join(np.insert(counts[i], 0, pgt).astype('string'))
                if not v_key in memo_hash:
                    v = log_multinomial(counts[i], lookup[pgt])
                    memo_hash[v_key] = v
                alt_list.append(memo_hash[v_key])

            memo_hash[key] = (r - max(alt_list))

        res[i] = memo_hash[key]
    return res


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