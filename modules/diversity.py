import numpy as np
import logging
import copy
import allel
import skbio
from pypopgen3.modules import treetools

logger = logging.getLogger()
logging.basicConfig(
    format='%(levelname)-8s %(asctime)s %(filename)  %(message)s')
#logging.basicConfig(format='%(levelname)-8s %(asctime)s %(funcName)20s()  %(message)s')
logger.setLevel(logging.WARNING)


def pairwise_differences(ac, an, correct_for_missing=True):
    """
    Calculate pairwise genetic difference between many
    individuals/populations/species from allele count matrix (2D-array).
    Allele counts can be from individual haplotypes (0 or 1),
    genotypes (0,1,2) (for diploids) or populations (0,1,2,..,an).
    an and ac have input dimensions number of variants x number of populations.
    If normalised by the accessible genome size, the diagonal of the output corresponds to within-population
    nucleotide diversity pi (or, for genotypes, heterozygosity), the off-diagonal entries correspond
    to nucleotide divergence d_xy.

    :param ac: 2D-array non-reference (or derived) allele count in haplotype/individual/population.
               Each row corresponds to a variant (e.g., SNP) site and each column corresponds to
               corresponds to a different haplotype/individual/population.
    :param an: 2D array of allele number for each variant in each haplotype/individual/population.
              Same dimensions as ac. Corresponds to the number of samples per haplotype/individual/population.
              E.g., an should be 1 for haplotypes and 2 for genotypes unless the variant call is missing.
              If the call is missing, then it is 0. For populations of diploid individuals an of SNP i in
              population j corresponds to the number of 2x the number of individuals in population j for
              which genotypes could be called at SNP i.
    :return: Symmetric (n populations x n populations) matrix of absolute counts of pairwise genetic difference.
             Each "population" can also just be a single individual or a single haplotype.

    """
    ac = np.array(ac)
    an = np.array(an).astype(float)
    an[an == 0] = np.nan
    af = ac / an
    af1 = np.nan_to_num(af)
    af1c = np.nan_to_num(1 - af)
    # this corrects for sampling without replacement in self comparisons
    af1c_within = np.nan_to_num(1 - (ac - 1.) / (an - 1))
    pi = np.zeros((af.shape[1], af.shape[1]))
    np.fill_diagonal(pi, (2 * af1 * af1c_within).sum(axis=0))
    dxy = np.einsum('ij,ik->jk', af1, af1c) + np.einsum('ij,ik->jk', af1c, af1)
    np.fill_diagonal(dxy, 0)
    pwd = pi + dxy

    if correct_for_missing:
        # should we use an for all or only for variable sites to correct for this?
        n_pwc = n_pairwise_comparisons(an)
        pwd = pwd * len(an) / n_pwc

    return pwd
#4

def n_pairwise_comparisons(an):
    has_data = np.array(an).astype(float)
    has_data[np.isnan(has_data)] = 0
    has_data[has_data > 0] = 1
    comparison_number = np.einsum('ij,ik->jk', has_data, has_data)
    return comparison_number

def pairwise_differences_from_gt_array(gt_arr, correct_for_missing=True):
    """
    Get pairwise genetic differences within and between genotypes from numpy array with genotpyes as loaded by
    allel.read_vcf.
    :param gt_arr: numpy array with genotpyes as loaded by allel.read_vcf
    :return: Symmetric (n individuals x n indiviudals) matrix of absolute counts of pairwise genetic difference.
             The diagonal corresponds to within-indiviudal heterozygosity.
    """
    gt = gt_arr.astype(float)
    gt[gt < 0] = np.nan
    gt = np.sum(gt, axis=-1)
    an = 2 * (gt == gt).astype(int)
    return pairwise_differences(gt, an, correct_for_missing=correct_for_missing)

def get_nj_tree_newick(pairwise_differences, samples):
    pairwise_differences = copy.deepcopy(pairwise_differences)
    np.fill_diagonal(pairwise_differences, 0)

    distance_matrix = skbio.DistanceMatrix(pairwise_differences,
                                           ids=samples)
    # This computes a neighbour-joining tree and outputs a newick string
    nj_newick = skbio.nj(distance_matrix, result_constructor=str)
    return nj_newick

def get_nj_tree(pairwise_differences, samples, outgroup=None, prune_outgroup=True):
    nj_newick = get_nj_tree_newick(pairwise_differences, samples)
    tree_inferred = treetools.HsTree(nj_newick)
    if outgroup is not None:
        tree_inferred.set_outgroup(outgroup)
        if prune_outgroup:
            tree_inferred.prune([n for n in tree_inferred.get_leaf_names() if n != outgroup])
    return tree_inferred

