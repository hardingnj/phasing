__author__ = 'Nicholas Harding'

from itertools import izip
from os import path, mkdir
import numpy as np


def get_fam_info(sampleid, ped_tbl):
    try:
        cross, role = ped_tbl.loc[sampleid]
    except KeyError:
        cross, role = 'none', 'none'
    if role != 'progeny':
        return cross, '0', '0'
    else:
        father, mother = ped_tbl[
            (ped_tbl.cross == cross) & (ped_tbl.role == 'parent')].index
        return cross, father, mother


def write_merlin(output_stem, genotypes, positions, samples, ped_tbl,
                 cm_per_base=2.75e-06, chrom='3L'):

    # ped file is horizontal... fam/id/mother father 0 0
    genotypes += 1
    ped_fh = open(output_stem + '.ped', 'w')

    for i, s in enumerate(samples):
        cross, father, mother = get_fam_info(s, ped_tbl)
        line = "\t".join([cross, s, father, mother, '0', '0'] +
                         np.apply_along_axis("/".join, 1,
                                             genotypes[:, i].astype(
                                                 'string')).tolist())
        ped_fh.write(line + "\n")
    ped_fh.close()

    # dat file
    dat_fh = open(output_stem + '.dat', 'w')
    dat_fh.write("\t".join(['T', 'phenotype']) + "\n")
    map_fh = open(output_stem + '.map', 'w')
    map_fh.write("\t".join(['CHROMOSOME', 'MARKER', 'POSITION']) + "\n")
    last_gp = 'X'

    for pos in positions:
        marker = 'id'+str(pos)
        dat_fh.write("\t".join(['M', marker]) + "\n")

        gp = str(round(cm_per_base*pos, 6))
        try:
            assert gp != last_gp
        except AssertionError:
            print(pos, last_gp)
            Exception(AssertionError)
        map_fh.write("\t".join([chrom, marker, gp]) + "\n")
    map_fh.close()
    dat_fh.close()


def write_duohmm(output_stem, genotypes, positions, samples, ped_tbl,
                 reference, alternate, chrom='3L'):

    if not path.isdir(output_stem):
        mkdir(output_stem)

    fh = open(output_stem + '.haps', 'w')
    for gt, pos, ref, alt in izip(genotypes, positions, reference, alternate):
        line = ' '.join([chrom, '.', str(pos), ref, alt] +
                        [str(q) for q in genotypes.flatten().tolist()])
        fh.write(line + "\n")
    fh.close()
    make_sample_file(samples, ped_tbl, output_stem + '.sample')


def make_sample_file(sample_list, ped_tbl, filepath):

    fh = open(filepath, 'w')
    fh.write('ID_1 ID_2 missing father mother sex plink_pheno' + "\n")
    fh.write('0 0 0 D D D B' + "\n")
    ped_counter = {p:0 for p in set(ped_tbl.cross)}

    for i, s in enumerate(sample_list):
        ped, father, mother = get_fam_info(s, ped_tbl=ped_tbl)
        ped_counter[ped] += 1
        fh.write(" ".join([ped + '_' + str(ped_counter[ped]), s, '0',
                           father, mother, '0', '-9']) + "\n")
    fh.close()