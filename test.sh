#!/usr/bin/env bash

# ----------------------
# Copyright 2023 PMG Lab
# Author: Lin Miao
# Licence: MIT
# Version: 20230825
# ----------------------

set -e

# Set the path of binaries
plink=plink
kggsee='java -Xmx20g -jar /app/pmglab/kggsee/kggsee.jar'
hess=hess.py
ldak=ldak


printf "\n\n1. Concatenate VCF file and convert to PLINK binary format ...\n\n"
cat metadata/auto_coding_with10kb_flk.vcf.gz.* >metadata/auto_coding_with10kb_flk.vcf.gz
rm metadata/auto_coding_with10kb_flk.vcf.gz.*

$plink --vcf metadata/auto_coding_with10kb_flk.vcf.gz --make-bed --out metadata/auto_coding_with10kb_flk
for chrom in {1..22}; do
  $plink --bfile metadata/auto_coding_with10kb_flk --make-bed --chr $chrom --out metadata/auto_coding_with10kb_flk.$chrom
done


printf "\n\n2a.Format test/Height.minimized.tsv.gz to the input formats of the three methods ...\n\n"
./format_sumstat.py --from-neale test/Height.minimized.tsv.gz --out-prefix test/Height

printf "\n\n2b.Format test/College.minimized.tsv.gz to the input formats of the three methods ...\n\n"
./format_sumstat.py --from-neale test/College.minimized.tsv.gz --out-prefix test/College

printf "\n\n2c.Format metadata/auto_coding_with10kb_flk.region.kggsee to the input formats of LDAK and HESS ...\n\n"
./format_region.py


printf "\n\n3a.Estimate h2 of Height by EHE ...\n\n"
$kggsee --gene-herit --filter-maf-le 0.01 --vcf-ref metadata/auto_coding_with10kb_flk.vcf.gz --regions-bed metadata/auto_coding_with10kb_flk.region.kggsee --sum-file test/Height.kggsee.gz --nmiss-col N --out test/Height.ehe

printf "\n\n3b.Estimate h2 of College by EHE ...\n\n"
$kggsee --gene-herit --filter-maf-le 0.01 --vcf-ref metadata/auto_coding_with10kb_flk.vcf.gz --regions-bed metadata/auto_coding_with10kb_flk.region.kggsee --sum-file test/College.kggsee.gz --nmiss-col N --out test/College.ehe


printf "\n\n4a.Estimate h2 of Height by HESS ...\n\n"
for chrom in {1..22}; do
  $hess --chrom $chrom --bfile metadata/auto_coding_with10kb_flk.$chrom --min-maf 0.01 --partition metadata/auto_coding_with10kb_flk.region.hess.$chrom --local-hsqg test/Height.hess.gz --out test/Height.hess.step1
done
$hess --tot-hsqg 0.5038 0.02279 --prefix test/Height.hess.step1 --out test/Height.hess.step2

printf "\n\n4b.Estimate h2 of College by HESS ...\n\n"
for chrom in {1..22}; do
  $hess --chrom $chrom --bfile metadata/auto_coding_with10kb_flk.$chrom --min-maf 0.01 --partition metadata/auto_coding_with10kb_flk.region.hess.$chrom --local-hsqg test/College.hess.gz --out test/College.hess.step1
done
$hess --tot-hsqg 0.2523 0.005995 --prefix test/College.hess.step1 --out test/College.hess.step2


printf "\n\n5a.Estimate h2 of Height by GBAT ...\n\n"
$ldak --cut-genes test/Height.gbat --genefile metadata/auto_coding_with10kb_flk.region.ldak --bfile metadata/auto_coding_with10kb_flk
$ldak --calc-genes-reml test/Height.gbat --summary test/Height.ldak --allow-ambiguous YES --bfile metadata/auto_coding_with10kb_flk --ignore-weights YES --power -0.25 &>test/Height.gbat/calc_genes_reml.log
$ldak --join-genes-reml test/Height.gbat

printf "\n\n5b.Estimate h2 of College by GBAT ...\n\n"
$ldak --cut-genes test/College.gbat --genefile metadata/auto_coding_with10kb_flk.region.ldak --bfile metadata/auto_coding_with10kb_flk
$ldak --calc-genes-reml test/College.gbat --summary test/College.ldak --allow-ambiguous YES --bfile metadata/auto_coding_with10kb_flk --ignore-weights YES --power -0.25 &>test/College.gbat/calc_genes_reml.log
$ldak --join-genes-reml test/College.gbat


printf "\n\n6.Collect all h2 estimates ...\n\n"
./collect_estimate.py --result-dir test --phenotypes Height,College --out test.result.tsv
