import itertools, pysam
from . import vcfpandas as vp


def get_within_between_comparisons(samples1, samples2, within_population_only_within_individual=True):
    '''
    within_population_only_within_individual... Within population, only use comparisons of the two haplotypes
                                                within individuals, not across indiviudals.
                                                This avoids phasing problems.

    Example:
        get_within_between_comparisons(['A','B'], ['C','D'])
    Returns:
        ('0-1,2-3',
         '4-5,6-7',
         '0-4,0-5,0-6,0-7,1-4,1-5,1-6,1-7,2-4,2-5,2-6,2-7,3-4,3-5,3-6,3-7')
    '''
    i = 0
    within1 = []
    within2 = []
    cross = []
    if within_population_only_within_individual:
        sep0 = '-'
    else:
        sep0 = ','
    for s in samples1:
        within1.append(str(2 * i) + sep0 + str(2 * i + 1))
        i += 1
    for s in samples2:
        within2.append(str(2 * i) + sep0 + str(2 * i + 1))
        i += 1

    for a, b in itertools.product(range(2 * len(samples1)),
                                  range(2 * len(samples1), 2 * (len(samples1) + len(samples2)))):
        cross.append(str(a) + '-' + str(b))
    return ','.join(within1), ','.join(within2), ','.join(cross)


def write_multihetsep(vcf_fn, outfile, accessible_bed, samples=None, chrom=None,
                      start=None, end=None, biallelic_only=True, keep_missing=True, debug=False):
    """
    Write a multihetsep file as used by msmc/msmc2.
    The tab separated file has the columns:
    chrom pos accessible_sites_since_last_positon haplotypes
    1	58432	63	TCCC
    1	58448	16	GAAA
    1	68306	15	CTTT
    See here for multihetsep file specification:
    https://github.com/stschiff/msmc/blob/master/guide.md

    Input:

    vcf_fn... a variant call format (vcf) file name. Usually containing biallelic SNPs.
              If indels are present they will be skipped.
    outfile... filename to write the hetsep file to.
    accessible_bed... .bed file with regions of the genome that could be called.
                            (non-filtered sites, this is the complement of a mask file)
    Optional arguments:
    samples... samples for which alleles/haplotypes should be written.
                (there will be two alleles per sample as this assumed diploid samples)
                Default: all samples in vcf
    chrom... chromosome which will be accessed (vcf required bcftools/tabix index)
    start... start position on chromosome (vcf required bcftools/tabix index)
    end... start position on chromosome (vcf required bcftools/tabix index)
    biallelic_only... skip sites with more than one alt allele
    keep_missing... keep sites where some indiviudals have missing genotypes
                    (missing genotypes are replaced by '?'
    debug... if True, only consider the first 1000 lines of the vcf (-region)
    """
    vcf_in = pysam.VariantFile(vcf_fn)
    last_pos = (start if start is not None else 0)

    if samples is None:
        samples = [s for s in vcf_in.header.samples]

    # open the output file for writing
    with open(outfile, 'w') as f:
        # iterate over lines/sites/SNPs in vcf
        for k, rec in enumerate(vcf_in.fetch(chrom, start, end)):
            # skip indels
            if len(rec.ref) > 1 or len(rec.alts[0]) > 1:
                continue
            # skip non biallelic
            if biallelic_only and len(rec.alts) > 1:
                continue

            # this will be a list of strings;
            # each string will contain two alleles per sample, e.g.: AATA
            # if het samples are unphased there will be one list entry per possible phasing,
            # e.g., ['AAAT','AATA']. If several samples are unphased there can be many combinations
            alleles = ['']

            # iterate over target samples
            for s in samples:
                # retrieve record for sample
                sr = rec.samples[s]
                # create tuple of sample alleles replacing missing calls by ?
                sample_alleles = tuple([(s.upper() if s is not None else '?') for s in sr.alleles])

                # skip if any missing allele
                if not keep_missing and '?' in sample_alleles:
                    continue

                # if sample is not phased and heterozygous
                if not sr.phased and sample_alleles[0] != sample_alleles[1]:
                    # make a copy of alleles
                    alleles1 = alleles[:]

                    # to each  string of alleles of previous samples
                    # add the alleles of current samples in reverse order
                    for i in range(len(alleles1)):
                        alleles1[i] += ''.join(sample_alleles[::-1])
                else:
                    # if phased, no additional orientations to add
                    alleles1 = []

                # to each  string of alleles of previous samples
                # add the alleles of current samples in phased order
                for i in range(len(alleles)):
                    alleles[i] += ''.join(sample_alleles)

                # add potentially additional allel combinations because of non-phased sample
                alleles = alleles + alleles1

            # if there is no variation across samples, skip this site
            if len(alleles) == 1 and len(set([a for a in alleles[0] if a != '?'])) == 1:
                continue

            # get number of accessible based between current and last site
            # accessible = max(1, vp.get_accessible_size(accessible_bed,
            #                                    rec.chrom, last_pos, rec.pos))
            accessible = vp.get_accessible_size(accessible_bed,
                                                rec.chrom, last_pos, rec.pos - 1) + 1

            # update last position to current
            last_pos = rec.pos
            # write the multihetsep line
            f.write('{}\t{}\t{}\t{}\n'.format(rec.chrom, rec.pos, accessible, ','.join(alleles)))
            if debug and k > 1000:
                break


