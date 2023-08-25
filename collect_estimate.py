#!/usr/bin/env python

# ----------------------
# Copyright 2023 PMG Lab
# Author: Lin Miao
# Licence: MIT
# Version: 20230825
# ----------------------

import re
import logging
import argparse
import numpy as np
import pandas as pd
from scipy.stats import norm

parser = argparse.ArgumentParser(description='Collect all heritability estimates into one file and convert the '
                                             'estimate of a dichotomous phenotype into liability-scale')
parser.add_argument('--result-dir', type=str, required=True, help='The result directory')
parser.add_argument('--phenotypes', type=str, required=True, help='Comma separated Abbr of phenotypes')
parser.add_argument('--out', type=str, required=True, help='The output file name')

args = parser.parse_args()
abbr_list = args.phenotypes.split(',')
logging.basicConfig(level=logging.DEBUG, format='[%(levelname)s] %(asctime)s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')


def cal_liabiliby_coef(k):
    if np.isnan(k):
        return 1
    else:
        return (k * (1 - k) / norm.pdf(norm.ppf(k))) ** 2 * 4


logging.info('Read metadata/phenotype.tsv and metadata/auto_coding_with10kb_flk.region.kggsee')
phenotype = pd.read_csv('metadata/phenotype.tsv', sep='\t', index_col='Abbr').loc[abbr_list, 'Prevalence']
liab_coef = phenotype.map(cal_liabiliby_coef)

gene_name = pd.read_csv('metadata/auto_coding_with10kb_flk.region.kggsee', sep='\t', header=None)
gene_name.set_index(3, drop=False, inplace=True)
gene_name[4] = np.int32((gene_name[1] + gene_name[2]) / 2)


logging.info('Collect EHE estimates')
df_dict = dict()
for abbr in abbr_list:
    df = pd.read_csv(f'{args.result_dir}/{abbr}.ehe.gene.pvalue.txt', sep='\t', index_col='Gene')[['Herit', 'HeritSE']]
    df_dict[abbr] = df

ehe_results = pd.concat(df_dict, axis=1)
ehe_results.index = pd.MultiIndex.from_frame(gene_name.loc[ehe_results.index, [3, 0, 4]])
ehe_results.index.names = [None, None, None]

ehe_results.loc[:, (slice(None), 'Herit')] = \
    (ehe_results.loc[:, (slice(None), 'Herit')].droplevel(1, axis=1) * liab_coef).values
ehe_results.loc[:, (slice(None), 'HeritSE')] = \
    (ehe_results.loc[:, (slice(None), 'HeritSE')].droplevel(1, axis=1) * liab_coef).values


logging.info('Collect HESS estimates')
df_dict = dict()
for abbr in abbr_list:
    df = pd.read_csv(f'{args.result_dir}/{abbr}.hess.step2.txt', sep='\t', index_col=['chr', 'start', 'end'])[
        ['local_h2g', 'se']]
    df_dict[abbr] = df

hess_results = pd.concat(df_dict, axis=1)
gene_name.set_index([0,1,2], drop=False, inplace=True)
hess_results.index = pd.MultiIndex.from_frame(gene_name.loc[hess_results.index, [3, 0, 4]])
hess_results.index.names = [None, None, None]

logging.info('Collect GBAT estimates')
reg = re.compile(r'Warning, Gene/Chunk (.+) is excluded')
df_dict = dict()
for abbr in abbr_list:
    exclude = list()
    with open(f'{args.result_dir}/{abbr}.gbat/calc_genes_reml.log') as I:
        for l in I:
            a = reg.match(l)
            if a:
                exclude.append(a.group(1))

    df = pd.read_csv(f'{args.result_dir}/{abbr}.gbat/remls.all', sep=' ', index_col='Gene_Name')[['Heritability', 'SD']]
    df[df.SD.isnull()] = np.nan
    df.loc[exclude] = np.nan
    df_dict[abbr] = df

gbat_results = pd.concat(df_dict, axis=1)
gbat_results.index = gbat_results.index.str.replace(r'\d+_', '', regex=True)
gene_name.set_index(3, drop=False, inplace=True)
gbat_results.index = pd.MultiIndex.from_frame(gene_name.loc[gbat_results.index, [3, 0, 4]])
gbat_results.index.names = [None, None, None]
gbat_results.loc[:, (slice(None), 'Heritability')] = \
    (gbat_results.loc[:, (slice(None), 'Heritability')].droplevel(1, axis=1) * liab_coef).values
gbat_results.loc[:, (slice(None), 'SD')] = \
    (gbat_results.loc[:, (slice(None), 'SD')].droplevel(1, axis=1) * liab_coef).values

logging.info(f'Write the results of all methods to {args.out}')
all_results = pd.concat({'EHE': ehe_results, 'HESS': hess_results, 'GBAT': gbat_results}, axis=1)
all_results = all_results.rename({'local_h2g': 'h2', 'Herit': 'h2', 'HeritSE': 'se', 'Heritability': 'h2', 'SD': 'se'},
                                 axis=1, level=2).swaplevel(0, 1, axis=1)[abbr_list]
all_results.index.names = ['HGNC_ID', 'CHR', 'MID_BP']
all_results.to_csv(args.out, sep='\t', na_rep='NA')
logging.info('Done.')
