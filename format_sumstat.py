#!/usr/bin/env python

# ----------------------
# Copyright 2023 PMG Lab
# Author: Lin Miao
# Licence: MIT
# Version: 20230825
# ----------------------

import logging
import argparse
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser(description='Format UKB GWAS results downloaded from the Neale Lab to '
                                             'the input formats of KGGSEE, HESS and LDAK')
parser.add_argument('--from-neale', type=str, required=True, help='GWAS result file downloaded from Neale Lab')
parser.add_argument('--out-prefix', type=str, required=True, help='Prefix of the outputs')
args = parser.parse_args()

logging.basicConfig(level=logging.DEBUG, format='[%(levelname)s] %(asctime)s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
logging.info('Read metadata/auto_coding_with10kb_flk.vcf.gz.')
bp_snp_map = np.loadtxt('metadata/auto_coding_with10kb_flk.vcf.gz', dtype=str, comments='#', usecols=[0, 1, 2, 3, 4])
bp_snp_map = pd.DataFrame(bp_snp_map[:, 2:], columns=['SNP', 'A2', 'A1'],
                          index=pd.MultiIndex.from_arrays(bp_snp_map[:, :2].T.astype(int)))

logging.info(f'Read {args.from_neale} and extract SNPs.')
from_neale = pd.read_csv(args.from_neale, sep='\t')
from_neale = from_neale[~from_neale.variant.str.match('X')]
from_neale.index = pd.MultiIndex.from_frame(from_neale.variant.str.split(':', expand=True).iloc[:, :2].astype(int))

df = pd.concat(
    [bp_snp_map, from_neale.loc[bp_snp_map.index, ['n_complete_samples', 'beta', 'se', 'tstat', 'pval']]], axis=1)
df.index.names = ['CHR', 'BP']
df.columns = ['SNP', 'A2', 'A1', 'N', 'BETA', 'SE', 'Z', 'P']

logging.info(f'Write a KGGSEE input file.')
df[['P', 'N']].to_csv(f'{args.out_prefix}.kggsee.gz', sep='\t')

logging.info(f'Write a HESS input file.')
df[['SNP', 'A1', 'A2', 'Z', 'N']].to_csv(f'{args.out_prefix}.hess.gz', sep='\t')

logging.info(f'Write an LDAK input file.')
df.rename({'N': 'n', 'SNP': 'Predictor'}, axis=1)[['Predictor', 'A1', 'A2', 'n', 'Z']]\
    .to_csv(f'{args.out_prefix}.ldak', sep='\t', index=False)

logging.info(f'Done.')
