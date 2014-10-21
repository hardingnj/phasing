#! /bin/bash

# simple bash script to filter the raw cross data by genotype qual
set -e
set -o pipefail

if test -z "$1"
	then echo "Error. GQ not specified"
	exit 1
fi

GQ=$1

echo $GQ
RAW=/data/anopheles/ag-crosses/data/release/0.1.GATK.PHASING.1000AG.AR2/h5/3L_ag-cross.h5
OUTDIR=/data/anopheles/ag1000g/data/1000g_09_13/evaluation/phasing/phase1.AR2/callsets/autosome/truth/vcf
BIN=/home/njh/phasing/bin/
PED=/home/njh/git/ag-crosses/meta/tbl_sample_ped.txt

OUT="${OUTDIR}/3L_ag-cross.${GQ}"

/home/njh/pyenv/science/bin/python ${BIN}/filter_hdf5.py ${RAW} ${OUT}.h5 -Q 30 -P ${PED}
/home/njh/pyenv/science/bin/python ${BIN}/hdf5_2_vcf.py ${OUT}.h5 ${OUT}.vcf
bgzip ${OUT}.vcf
tabix -p vcf ${OUT}.vcf.gz
touch ${OUT}.vcf.gz.ok
