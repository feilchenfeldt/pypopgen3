import numpy as np
import pandas as pd
import logging
import copy
import allel
import scipy, skbio, warnings
from pypopgen3.modules import treetools

#import warnings
#warnings.filterwarnings("ignore",category=UserWarning)

import logging
# logger = logging.getLogger()
# logging.basicConfig(
#     format='%(levelname)-8s %(asctime)s %(filename)  %(message)s')
# logging.basicConfig(format='%(levelname)-8s %(asctime)s %(funcName)20s()  %(message)s')
# logger.setLevel(logging.WARNING)
logging.getLogger("allel").setLevel(logging.ERROR)

def get_hap_pwd_arr(vcf_fn, chrom=None, start=None, end=None, samples=None, ancestral_sample=None):
    region = ''

    if chrom is not None:
        region += chrom
    if start is not None:
        assert chrom is not None
        region += ':' + str(int(start))
    if end is not None:
        assert start is not None
        region += '-' + str(int(end))
    if not region:
        region = None

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="invalid INFO header:")
        callset = allel.read_vcf(vcf_fn, fields=['calldata/GT'],  # samples=samples,
                             alt_number=1, region=region, tabix='tabix')
    if callset is not None:
        try:
            gt_arr = callset['calldata/GT']
        except:
            print(callset)
        # gt = allel.GenotypeArray(gt_arr)
        ht_arr = gt_arr.reshape(gt_arr.shape[0], -1)

        if (ht_arr == -1).any().any():
            if ancestral_sample is None:
                # for now, assume missing calls are ref, will be true most of the time
                # workaround "bug" in allel, where it uses the literal -1 value for missing genotypes when calculating distance
                ht_arr[ht_arr == -1] = 0
            else:
                # assume missing calls are ancestral
                _, _, _, _, samples = allel.read_vcf_headers(vcf_fn)
                ancestral_ix = samples.index(ancestral_sample)
                given_axis = 0
                h = ht_arr[:,ancestral_ix]
                # assume if ancestral state is missing that ref is ancestral
                h[h == -1] = 0
                ht_arr[:,ancestral_ix] = h
                a = (ht_arr == -1)
                b = ht_arr[:, ancestral_ix * 2]
                dim_array = np.ones((1, a.ndim), int).ravel()
                dim_array[given_axis] = -1
                b_reshaped = b.reshape(dim_array)
                missing_replaced_by_anc = a * b_reshaped
                ht_arr = ht_arr + (ht_arr == -1) + missing_replaced_by_anc

        ht = allel.HaplotypeArray(ht_arr)
        d = allel.pairwise_distance(ht, metric='cityblock')
    else:
        if samples is None:
            _, _, _, _, samples = allel.read_vcf_headers(vcf_fn)
        d = np.zeros(2 * len(samples) ** 2 - len(samples))
    return d


def get_samples(vcf_fn=None, samples=None):
    assert vcf_fn is not None or samples is not None, "at least one argument must be given"
    if samples is None:
        _, _, _, _, samples = allel.read_vcf_headers(vcf_fn)
    return samples


def hap_pwd_to_df(pwd_arr, samples):
    d_mat = scipy.spatial.distance.squareform(pwd_arr)
    ix = pd.MultiIndex.from_product([samples, [0, 1]])
    pwd_hap = pd.DataFrame(d_mat, index=ix, columns=ix)
    return pwd_hap


def hap_to_dip(pwd_hap_df):
    pwd_diplo_allel = pwd_hap_df.groupby(axis=1, level=0).sum().groupby(axis=0, level=0).sum() / 4
    pwd_diplo_allel = pwd_diplo_allel + np.diag(np.diag(pwd_diplo_allel))
    #Make symmetric, check why this is needed
    pwd_diplo_allel = (pwd_diplo_allel + pwd_diplo_allel.T) / 2
    return pwd_diplo_allel


def get_hap_pwd(vcf_fn, chrom=None, start=None, end=None, samples=None, ancestral_sample=None):
    hap_pwd_arr = get_hap_pwd_arr(vcf_fn, chrom, start, end, samples, ancestral_sample)
    samples = get_samples(vcf_fn=vcf_fn, samples=samples)
    hap_df = hap_pwd_to_df(hap_pwd_arr, samples)
    return hap_df


def get_dip_pwd(vcf_fn, chrom=None, start=None, end=None, samples=None, ancestral_sample=None):
    hap_df = get_hap_pwd(vcf_fn, chrom, start, end, samples, ancestral_sample)
    dip_df = hap_to_dip(hap_df)
    return dip_df


def get_total_pwd(vcf_fn, chrom=None, start=None, end=None, chunksize=50000, samples=None, ancestral_sample=None):
    hap_pwd_arr = get_hap_pwd_arr(vcf_fn, chrom, start, min(start + chunksize - 1, end), samples, ancestral_sample)

    for s in range(start + chunksize, end, chunksize):
        e = min(s + chunksize - 1, end)
        hap_pwd_arr += get_hap_pwd_arr(vcf_fn, chrom, s, e, samples, ancestral_sample)
    samples = get_samples(vcf_fn=vcf_fn, samples=samples)

    return hap_pwd_to_df(hap_pwd_arr, samples)


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
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="invalid value encountered in true_divide")
        warnings.filterwarnings("ignore", message="divide by zero encountered in true_divide")
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
#4dfs

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
        tree_inferred.set_outgroup(outgroup, end_at_present=False)
        if prune_outgroup:
            tree_inferred.prune([n for n in tree_inferred.get_leaf_names() if n != outgroup])
    return tree_inferred

