#! /bin/bash
# create the merlin files: ped/dat/map
set -e
set -o pipefail

# This is a bash script that does much of the heavy lifting for the creation of VCF files
PYENV=/home/njh/pyenv/science/bin/python
PYBIN=/home/njh/git/phasing/bin
PEDTB=/home/njh/git/ag-crosses/meta/tbl_sample_ped.txt
MINGQ=40

OUTDI=/data/anopheles/ag1000g/data/1000g_09_13/evaluation/phasing/phase1.AR2/callsets/autosome/truth

WORKD=${OUTDI}/gq${MINGQ}
FSTEM=3L_ag-cross.gq${MINGQ}

# find all vcfs
crosses=( 18-5 22-1 29-2 36-9 37-3 42-4 45-1 46-9 47-6 56-1 61-3 73-2 78-2
80-2 )

for cr in ${crosses[@]}
do $PYENV $PYBIN/create_merlin_files.py ${WORKD}/${FSTEM}_${cr}.vcf.gz \
  ${WORKD}/${FSTEM}_${cr}.sample ${WORKD}/merlin
done