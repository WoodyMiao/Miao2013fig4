#!/usr/bin/env python

# ----------------------
# Copyright 2023 PMG Lab
# Author: Lin Miao
# Licence: MIT
# Version: 20230825
# ----------------------

import pandas as pd

df = pd.read_csv('metadata/auto_coding_with10kb_flk.region.kggsee', sep='\t', header=None)
df[[3, 0, 1, 2]].to_csv('metadata/auto_coding_with10kb_flk.region.ldak', sep='\t', header=False, index=False)

df.columns = ['chr', 'start', 'stop', 'gene']
df['chr'] = 'chr' + df['chr'].astype(str)
for c in range(1, 23):
    df.loc[df['chr'] == f'chr{c}', ['chr', 'start', 'stop']].to_csv(f'metadata/auto_coding_with10kb_flk.region.hess.{c}', sep='\t', index=False)
