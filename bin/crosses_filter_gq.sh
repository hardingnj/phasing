#! /bin/bash

# simple bash script to filter the raw cross data by genotype qual
set -e
set -o pipefail

if test -z "$1"
	then echo "Error. GQ not specified"
	exit 1
fi

GQ=$1

RAW=/data/anopheles/ag-crosses/data/release/0.1.GATK.PHASING.1000AG.AR2/3L_ag-cross.vcf.gz
OUTDIR=/data/anopheles/ag1000g/data/1000g_09_13/evaluation/phasing/phase1.AR2/callsets/autosome/truth/vcf

OUT="${OUTDIR}/3L_ag-cross.${GQ}"

vcftools --gzvcf ${RAW} --minGQ ${GQ} --recode --out ${OUT}
bgzip ${OUT}.recode.vcf
tabix -p vcf ${OUT}.recode.vcf.gz
touch ${OUT}.recode.vcf.gz.ok



