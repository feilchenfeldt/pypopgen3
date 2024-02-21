import subprocess
import numpy as np
import pandas as pd
#import logging
import copy
import allel
import scipy, skbio, warnings
from pypopgen3.modules import treetools
from pypopgen3.modules import vcfpandas

#import warnings
#warnings.filterwarnings("ignore",category=UserWarning)

#import logging
# logger = logging.getLogger()
# logging.basicConfig(
#     format='%(levelname)-8s %(asctime)s %(filename)  %(message)s')
# logging.basicConfig(format='%(levelname)-8s %(asctime)s %(funcName)20s()  %(message)s')
# logger.setLevel(logging.WARNING)
#logging.getLogger("allel").setLevel(logging.ERROR)
#
#print('test')

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
        warnings.filterwarnings("ignore", message=".*invalid INFO header:.*")
        warnings.filterwarnings('error', '.*Could not retrieve index.*', )
        callset = allel.read_vcf(vcf_fn, fields=['calldata/GT'],   samples=samples,
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
                if samples is None:
                    _, _, _, _, samples = allel.read_vcf_headers(vcf_fn)
                ancestral_ix = samples.index(ancestral_sample)
                given_axis = 0
                try:
                    h = ht_arr[:,ancestral_ix]
                except Exception as e:
                    print(len(samples),ht_arr.shape)
                    raise e
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


def get_samples(vcf_fn, samples=None):
    #assert vcf_fn is not None or samples is not None, "at least one argument must be given"
    _, _, _, _, vcf_samples = allel.read_vcf_headers(vcf_fn)
    if samples is not None:
        samples1 = [s for s in samples if s in vcf_samples]
    else:
        samples1 = vcf_samples

    return samples1


def hap_pwd_to_df(pwd_arr, samples):
    d_mat = scipy.spatial.distance.squareform(pwd_arr)
    ix = pd.MultiIndex.from_product([samples, [0, 1]])
    pwd_hap = pd.DataFrame(d_mat, index=ix, columns=ix)
    return pwd_hap


def hap_to_dip(pwd_hap_df):
    try:
        pwd_diplo_allel = pwd_hap_df.groupby(axis=1, level=0).sum().groupby(axis=0, level=0).sum() / 4
    except TypeError:
        print('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX')
        dip = pwd_hap_df.groupby(axis=1, level=0).sum().groupby(axis=0, level=0).sum()
        print(dip)
        print(dip.stack().apply(type).value_counts())
        print('ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ')
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
    samples = get_samples(vcf_fn, samples)
    dip_df = hap_to_dip(hap_df)
    return dip_df.loc[samples, samples]

def get_local_tree(vcf_fn, chrom, start, end, samples=None, ancestral_sample=None):
    pwd =  get_dip_pwd(vcf_fn, chrom=chrom, start=start, end=end, samples=samples, ancestral_sample=ancestral_sample)
    samples = get_samples(vcf_fn, samples)
    tree = get_nj_tree(pwd.loc[samples, samples].values, samples, outgroup=ancestral_sample, prune_outgroup=True)
    return tree


def get_total_pwd(vcf_fn, chrom, start, end, chunksize=50000, samples=None, ancestral_sample=None):

    assert start is not None and end is not None, "start/end must be explicitly given"

    end0 = min(start + chunksize - 1, end)
    hap_pwd_arr = get_hap_pwd_arr(vcf_fn, chrom, start, end0, samples, ancestral_sample)

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

def get_nj_tree(pairwise_differences, samples, outgroup=None,outgroups=None, prune_outgroup=True):


    pairwise_differences = np.array(pairwise_differences)
    np.fill_diagonal(pairwise_differences, 0)

    distance_matrix = skbio.DistanceMatrix(pairwise_differences,
                                           ids=samples)
    # This computes a neighbour-joining tree and outputs a newick string
    tree_inferred = skbio.nj(distance_matrix, result_constructor=treetools.HsTree)

    assert outgroup is None or outgroups is None, "cannot provide both outgroup and outgroups"
    if outgroup is not None or outgroups is not None:
        if outgroups is not None:
            outgroup = outgroups[0]
        elif outgroup is not None:
            outgroups = [outgroup]
        tree_inferred.set_outgroup(outgroup, end_at_present=False)
        if prune_outgroup:
            tree_inferred.prune([n for n in tree_inferred.get_leaf_names() if n not in outgroups])
    return tree_inferred


def get_pwd_and_trees(vcf_fn, chrom, start, end, samples=None, outgroup=None):
    """
    Calculates pairwise differences and NJ tree for a genomic region.

    Returns:
    haplotype tree (two tips per samples).
    sample tree
    sample pairwise difference matrix (as pandas data frame)
    """

    samples = get_samples(vcf_fn, samples)
    pwd_arr = get_hap_pwd_arr(vcf_fn, chrom=chrom, start=start, end=end,
                                 samples=samples, ancestral_sample=outgroup)
    hap_pwd_df = hap_pwd_to_df(pwd_arr, samples).fillna(0)

    hap_samples = ["_".join([t[0], str(t[1])]) for t in hap_pwd_df.index]
    if outgroup is not None:
        outgroups = [outgroup + '_0', outgroup + '_1']
    else:
        outgroups = None

    hap_tree = get_nj_tree(hap_pwd_df.values, hap_samples,
                              outgroup=None, outgroups=outgroups, prune_outgroup=True)
    dip_pwd_df = hap_to_dip(hap_pwd_df)

    dip_tree = get_nj_tree(dip_pwd_df.values, dip_pwd_df.index,
                              outgroup=outgroup, outgroups=None, prune_outgroup=True)

    return hap_tree.write()[0], dip_tree.write()[0], dip_pwd_df


def get_ac_an(gt_str):
    try:
        return int(gt_str[0]) +int(gt_str[2]), 2
    except ValueError:
        return 0, 0

def nj_tree_from_vcf(vcf_fn, chrom, start, end, samples=None, outgroup=None, prune_outgroup=True):
    if samples is None:
        header, _ = vcfpandas.parse_vcf_header(vcf_fn)
        samples = header[9:]
    sample_str = ','.join(samples)
    vcf_fn0 = vcf_fn#.format(chrom)
    c = f"bcftools view -H --samples {sample_str} --regions {chrom}:{start}-{end} {vcf_fn0}"
    p = subprocess.Popen(c, shell=True, stdout=subprocess.PIPE)
    #bcftools query -l input.vcf
    ac_an = []
    for line in p.stdout.readlines():
        ac_an.append([get_ac_an(s) for s in line.decode().split('\t')[9:]])
    ac_an = np.array(ac_an)
    ac = ac_an[:,:,0]
    an = ac_an[:,:,1]
    pwd = pairwise_differences(ac,an)
    if np.any(np.isnan(pwd)):
        print('Warning: nans in distance matrix. Setting to Zero')
        pwd[np.isnan(pwd)] = 0


    tree = get_nj_tree(pwd, samples, outgroup=outgroup, prune_outgroup=prune_outgroup)

    return tree


def is_supported(node, test_tree):
    leaves = node.get_leaf_names()
    b_node = test_tree.get_common_ancestor(leaves)
    b_leaves = b_node.get_leaf_names()
    if set(leaves) != set(b_leaves):
        return False
    else:
        return True


def get_bootstrap_support(real_tree, pwd_window, n_bootstrap_samples):
    support_dic = {}

    support_tree = copy.deepcopy(real_tree)
    n_windows = pwd_window.shape[1]

    for node in support_tree.iter_descendants():
        if not node.is_leaf():
            setattr(node, 'n_support', 0)
    for i in range(n_bootstrap_samples):
        t = np.random.choice(n_windows, size=n_windows)
        sample = pwd_window.iloc[:, t]
        sample_pwd = sample.sum(axis=1).unstack()

        # convert the pairwise differences to an object used by skbio
        distance_matrix = skbio.DistanceMatrix(sample_pwd,
                                               ids=sample_pwd.index)
        # This computes a neighbour-joining tree and outpus a newick string
        nj_newick = skbio.nj(distance_matrix, result_constructor=str)
        # load the tree into a tree object
        # (this object class is defined in Hannes Svardal packaged
        # pypogen3. It is very similar to an ete3 tree object, but
        # has a useful plot method implemented.
        test_tree = treetools.HsTree(nj_newick)
        test_tree.set_outgroup(outgroup, end_at_present=False)

        for node in support_tree.iter_descendants():
            if not node.is_leaf():
                supported = is_supported(node, test_tree)
                node.n_support += int(supported)
    for node in support_tree.iter_descendants():
        if not node.is_leaf():
            setattr(node, 'pct_support', 100 * node.n_support / n_bootstrap_samples)
            support_dic.update({node.get_name(): node.pct_support})
    support_s = pd.Series(support_dic)
    support_s.name = "Percentage bootstrap support"
    return support_tree, support_s