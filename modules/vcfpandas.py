"""
Tools that help to load/parse information
from a VCF (variant call format) file into 
pandas data frames and series.
Tools to handle such data.
"""

import re
import gc
import subprocess
import itertools
import logging
import gzip
import numpy as np
import pandas as pd

logger = logging.getLogger()
#logging.basicConfig(
#    format='%(levelname)-8s %(asctime)s %(filename)  %(message)s')
#logging.basicConfig(format='%(levelname)-8s %(asctime)s %(funcName)20s()  %(message)s')
#logger.setLevel(logging.WARNING)

def get_accessible_size(accessible_fn, chrom, start=None, end=None):
    """
    Returns accessible size in bp of a given genomic region from a tabix
    indexed bed.gz file.

    Requires tabix to be installed and in the PATH.
    """
    region_str = str(chrom)
    if start is not None:
        assert end is not None, "start and end most both be given or none"
        region_str += f":{start}-{end}"
    elif end is not None:
        raise Exception("start and end most both be given or none")

    p = subprocess.Popen(f"tabix {accessible_fn} {region_str}", shell=True, stdout=subprocess.PIPE)
    d = pd.read_csv(p.stdout, sep='\t', header=None,
                    names=['chrom', 'start', 'end'], usecols=['start', 'end'])
    # check that intervals are non overlapping
    assert (d.iloc[:-1]['end'].values <= d.iloc[1:]['start'].values).all()
    try:
        d.iloc[0]['start'] = start
        d.iloc[-1]['end'] = end
        access = (d['end'] - d['start']).sum()
    except IndexError:
        access = 0
    return access


def add_info_line(info_dic, line):
    """
    This function parses a vcf header line
    and adds the information in the line
    into info_dic. 
    """
    # split all commas except if in quotes
    r = re.compile('[^,]+".+"[^,]*'
                   '|'
                   '[^,]+')
    try:
        category, value = line.strip().split('=<')
    except ValueError:
        try:
            info_dic["uncategorised"].update(
                {line.strip().split('=')[0][2:]: line.strip().split('=')[1]})
        except KeyError:
            try:
                info_dic.update({"uncategorised": {line.strip().split('=')[
                                0][2:]: line.strip().split('=')[1]}})
            except Exception as e:
                print(line)
                print(str(e))
                logging.warning("Cannot parse {}".format(line))
                logging.warning(str(e))
                try:
                    info_dic["unparsable"].update(line.strip())
                except KeyError:
                    info_dic.update({"unparsable":line.strip()})
        except IndexError:
            try:
                info_dic["unparsable"].update(line.strip())
            except KeyError:
                info_dic.update({"unparsable":line.strip()})
        return
    category = category[2:]
    tags = r.findall(value[:-1])
    subdic = {k: v for k, v in [t.split('=', 1) for t in tags]}
    try:
        info_dic[category].update({subdic["ID"]: subdic})
    except KeyError:
        info_dic.update({category: {subdic["ID"]: subdic}})


def parse_vcf_header(fn):
    """
    This function parses the header of a vcf file
    and returns a list of the table header (last header line)
    and a dictionary of the other header entries.
    """
    info_dic = {}
    fh = gzip.open(fn) if fn[-3:] == ".gz" else open(fn)
    for line in fh:
        line = line.decode('utf-8')
        if line[:6] == "#CHROM":
            header = line.strip().split('\t')
            header[0] = header[0][1:]
        else:
            add_info_line(info_dic, line)
        if line[0] != '#':
            break
    return header, info_dic


def map_reduce_haplo(fn, fun, samples_h0=None, samples_h1=None, chrom=None,
                     start=None, end=None, chunksize=50000, apply_fun=None,
                     get_result_fun=None, reduce_fun=None, **read_args):
    """
    Read haplotypes from specific region tabixed vcf.gz file
    and calculate pairwise differences.

    """
    if apply_fun is None:
        def apply_fun(fun, chunk0, chunk1):
            return fun(chunk0, chunk1)
    if reduce_fun is None:
        reduce_fun = lambda res: res

    header, _ = parse_vcf_header(fn)

    if samples_h0 is None and samples_h1 is None:
        samples_h0 = header[9:]
        samples_h1 = header[9:]
    elif samples_h0 is None:
        samples_h0 = []
        samples_h1 = [str(s) for s in samples_h1]
    elif samples_h1 is None:
        samples_h1 = []
        samples_h0 = [str(s) for s in samples_h0]

    t0 = get_vcf_df(fn, chrom=chrom, start=start, end=end, header=header, usecols=['CHROM', 'POS'] + samples_h0,
                    converters=converters.first_haplotype(samples_h0), chunksize=chunksize, **read_args)
    t1 = get_vcf_df(fn, chrom=chrom, start=start, end=end, header=header, usecols=['CHROM', 'POS'] + samples_h1,
                    converters=converters.second_haplotype(samples_h1), chunksize=chunksize, **read_args)

    results = []
    for chunk0, chunk1 in itertools.izip(t0, t1):
        results.append(apply_fun(fun, chunk0, chunk1))
        gc.collect()

    # reduce
    result = reduce_fun([get_result_fun(r) for r in results]
                        if get_result_fun is not None else results)

    return result


def map_fly_reduce_haplo(fn, fun, samples_h0=None, samples_h1=None, chrom=None,
                     start=None, end=None, chunksize=50000, apply_fun=None,
                     get_result_fun=None, fly_reduce_fun=None, **read_args):
    """
    Read haplotypes from specific region tabixed vcf.gz file
    and apply function.

    Reduces on the fly by applying reduce function after each 
    chunk.

    """
    if apply_fun is None:
        def apply_fun(fun, chunk0, chunk1):
            return fun(chunk0, chunk1)
    if fly_reduce_fun is None:
        def fly_reduce_fun(chunk_res, result=None):
            if result is None:
                return chunk_res
            else:
                return result + chunk_res

    header, _ = parse_vcf_header(fn)

    if samples_h0 is None and samples_h1 is None:
        samples_h0 = header[9:]
        samples_h1 = header[9:]
    elif samples_h0 is None:
        samples_h0 = []
        samples_h1 = [str(s) for s in samples_h1]
    elif samples_h1 is None:
        samples_h1 = []
        samples_h0 = [str(s) for s in samples_h0]
    #try: 
    t0 = get_vcf_df(fn, chrom=chrom, start=start, end=end, header=header, usecols=['CHROM', 'POS'] + samples_h0,
                    converters=converters.first_haplotype(samples_h0), chunksize=chunksize, **read_args)
    #except ValueError:
    #    print header
    #    print ['CHROM', 'POS'] + samples_h0

    t1 = get_vcf_df(fn, chrom=chrom, start=start, end=end, header=header, usecols=['CHROM', 'POS'] + samples_h1,
                    converters=converters.second_haplotype(samples_h1), chunksize=chunksize, **read_args)

    result = None
    for chunk0, chunk1 in itertools.izip(t0, t1):
        chunk0.columns = pd.MultiIndex.from_arrays(
                [chunk0.columns, [0] * len(chunk0.columns)])
        chunk1.columns = pd.MultiIndex.from_arrays(
                [chunk1.columns, [1] * len(chunk1.columns)])
        result = fly_reduce_fun(apply_fun(fun, chunk0, chunk1), result)
        gc.collect()

    return result


def map_reduce_geno(fn, fun, samples=None, chrom=None, start=None, end=None,
                    chunksize=50000, apply_fun=None, get_result_fun=None,
                    reduce_fun=None, **read_args):
    """

    """
    if apply_fun is None:
        apply_fun = lambda fun, chunk: fun(chunk)
    if reduce_fun is None:
        reduce_fun = lambda res: res

    header, _ = parse_vcf_header(fn)

    if samples is None:
        samples = header[9:]
    else:
        samples = [str(s) for s in samples]

    t = get_vcf_df(fn, chrom=chrom, start=start, end=end, header=header, usecols=['CHROM', 'POS'] + samples,
                   converters=converters.genotype_converter(samples), chunksize=chunksize, **read_args)

    results = []
    for chunk in t:
        results.append(apply_fun(fun, chunk))
        gc.collect()

    # reduce
    result = reduce_fun([get_result_fun(r) for r in results]
                        if get_result_fun is not None else results)

    return result


def get_vcf_df(fn, chrom=None, start=None, end=None, header=None, **read_csv_args):
    """
    This function reads a region from a bgzipped and 
    tabix indexed VCF into a pandas data frame.
    """

    if header is None:
        header, _ = parse_vcf_header(fn)

    if fn[-3:] == ".gz" and chrom is not None:
        region = str(chrom)
        if end is not None and start is None:
            start = 0
        if start is not None:
            region += ':' + str(start)
            if end is not None:
                region += '-' + str(end)

        tabix_stream = subprocess.Popen(['tabix', fn, region],
                                        stdout=subprocess.PIPE,
                                        stderr=subprocess.PIPE)
        vcf_df = pd.read_csv(tabix_stream.stdout, sep="\t",
                             index_col=[0, 1], header=None,
                             names=header, **read_csv_args)

    else:
        assert chrom is None, ("Custom chromosome only supported for bgzipped and"
                               "tabixed .gz files")
        assert start is None, ("Custom start only supported for bgzipped and"
                               "tabixed .gz files")
        assert end is None, ("Custom end only supported for bgzipped and"
                             "tabixed .gz files")
        vcf_df = pd.read_csv(fn, sep="\t",
                             index_col=[0, 1], header=None,
                             names=header, comment='#', **read_csv_args)

    return vcf_df


def get_genotype(s):
    """
    Get genotype from vcf genotype string.
    Attention: Only parses biallelic genotypes 
    at the moment.
    """
    gt_str = s.split(':', 1)[0]
    if gt_str[:3] in ["0/0", "0|0"]:
        return 0
    elif gt_str[:3] in ["1/1", "1|1"]:
        return 2
    elif gt_str[:3] in ["0/1", "0|1", "1|0"]:
        return 1
    else:
        return np.nan


def get_first_haplotype(s):
    gt_str = s.split(':', 1)[0]
    try:
        return int(gt_str[0])
    except ValueError:
        assert gt_str[0] == '.', 'Unknown allele {}'.format(gt_str[0])
        return np.nan
    #except IndexError:
        #this should at least lead to a clear warning. It points to 
        #an error in the vcf formatting!!!!
    #    return np.nan


def get_second_haplotype(s):
    gt_str = s.split(':', 1)[0]
    try:
        return int(gt_str[2])
    except ValueError:
        assert gt_str[2] == '.', 'Unknown allele {}'.format(gt_str[2])
        return np.nan
    #except IndexError:
        #print gt_str
    #    return np.nan


def get_depth(s):
    try:
        depth_str = s.split(':')[2]
    except IndexError:
        return np.nan
    try:
        return int(depth_str)
    except ValueError:
    #    assert gt_str[2] == '.', 'Unknown allele {}'.format(gt_str[2])
        return np.nan


def get_haplotype(s):
    gt_str = s.split(':', 1)[0]
    try:
        return int(gt_str[0]), int(gt_str[2])
    except ValueError:
        assert gt_str[0] == '.', 'Unknown allele {}'.format(gt_str[0])
        assert gt_str[2] == '.', 'Unknown allele {}'.format(gt_str[2])
        return np.nan, np.nan


def get_info_dic(s):
    """
    Parse the string from the VCF INFO
    column and return it as a dictionary.
    """
    info_tuples = [t.split('=') for t in s.split(';')]
    info_tuples = [t for t in info_tuples if len(t) == 2]
    tags = [t for t in info_tuples if len(t) != 2]
#    try:
    d = {k: v for (k, v) in info_tuples}
    d.update({'_tags': tags})
#    except ValueError:
#        d = {}
#        logging.warning("Can not parse info column at {}:{}".format(line[0],line[1]))
    return d


class converters:
    """
    This is just a container for 
    different static converter methods.
    Not meant to be instantiated.
    """
    @staticmethod
    def genotype_converter(samples):
        return {n: get_genotype for n in samples}

    @staticmethod
    def haplotype_converter(samples):
        return {n: get_haplotype for n in samples}

    @staticmethod
    def first_haplotype(samples):
        return {n: get_first_haplotype for n in samples}

    @staticmethod
    def second_haplotype(samples):
        return {n: get_second_haplotype for n in samples}

    @staticmethod
    def depth(samples):
        return {n: get_depth for n in samples}

def get_genotype_df(fn, chrom=None, start=None, end=None, header=None, samples=None):

    if header is None:
        header, _ = parse_vcf_header(fn)


    if samples is None:
        samples = header[9:]
    else:
        samples = [str(s) for s in samples]

    vcf_df = get_vcf_df(fn, chrom, start, end,
                        header=header,
                        usecols=["CHROM", "POS"] + samples,
                        converters=converters.genotype_converter(samples))

    return vcf_df


def get_haplotype_df(fn,  chrom=None,
                     start=None, end=None, samples=None, samples_h0=None, samples_h1=None,  **read_args):
    """
    Read haplotypes from specific region tabixed vcf.gz file.

    """
    header, _ = parse_vcf_header(fn)

    if samples is not None:
        assert (samples_h0 is None) and (samples_h1 is None)
        samples_h0 = samples
        samples_h1 = samples
    elif samples_h0 is not None or samples_h1 is not None:
        samples_h0 = samples_h0 if samples_h0 is not None else []
        samples_h1 = samples_h1 if samples_h1 is not None else []
    else:
        samples_h0 = header[9:]
        samples_h1 = header[9:]
    #    assert (samples_h0 is None) == (samples_h1 is None)
    #    if samples_h1 is None:        
    #        samples_h0 = header[9:]
    #        samples_h1 = header[9:]


    t0 = get_vcf_df(fn, chrom=chrom, start=start, end=end, header=header, usecols=['CHROM', 'POS'] + samples_h0,
                    converters=converters.first_haplotype(samples_h0),  **read_args)
    t1 = get_vcf_df(fn, chrom=chrom, start=start, end=end, header=header, usecols=['CHROM', 'POS'] + samples_h1,
                    converters=converters.second_haplotype(samples_h1), **read_args)

    # t0.columns = pd.MultiIndex.from_arrays(
    #         [t0.columns, [0] * len(t0.columns)])
    # t1.columns = pd.MultiIndex.from_arrays(
    #         [t1.columns, [1] * len(t1.columns)])
    # hap_df  = pd.concat([t0, t1], axis=1).sort_index(axis=1)

    t0.columns = pd.MultiIndex.from_arrays(
            [t0.columns, [0] * len(t0.columns)])
    t1.columns = pd.MultiIndex.from_arrays(
            [t1.columns, [1] * len(t1.columns)])

    hap_df  = pd.concat([t0, t1], axis=1).sort_index(axis=1)
    hap_df = hap_df[ samples_h0 + [s for s in samples_h1 if s not in samples_h0]]


    return hap_df


def get_accessible_size(accessible_fn, chrom, start=None, end=None):
    """
    Returns accessible size in bp of a given genomic region from a tabix
    indexed bed.gz file.


    """
    region_str = str(chrom)
    if start is not None:
        assert end is not None, "start and end most both be given or none"
        region_str += f":{start}-{end}"
    elif end is not None:
        raise Exception("start and end most both be given or none")


    p = subprocess.Popen(f"tabix {accessible_fn} {region_str}", shell=True, stdout=subprocess.PIPE)
    d = pd.read_csv(p.stdout, sep='\t', header=None,
                    names=['chrom', 'start', 'end'], usecols=['start', 'end'])
    # check that intervals are non overlapping
    assert (d.iloc[:-1]['end'].values <= d.iloc[1:]['start'].values).all()
    try:
        if start is not None:
            d.iloc[0]['start'] = max(d.iloc[0]['start'], start)
        if end is not None:
            d.iloc[-1]['end'] = min(d.iloc[-1]['end'], end)
        access = (d['end'] - d['start']).sum()
    except IndexError:
        access = 0
    return access
