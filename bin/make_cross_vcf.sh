#! /bin/bash
set -e
set -o pipefail

# This is a bash script that does much of the heavy lifting for the creation of VCF files
PYENV=/home/njh/pyenv/science/bin/python
PYBIN=/home/njh/git/phasing/bin
VKEEP=/home/njh/exec/bin/vcfkeepsamples
RAWH5=/data/anopheles/ag-crosses/data/release/0.1.GATK.PHASING.1000AG.AR2/h5/3L_ag-cross.h5 
PEDTB=/home/njh/git/ag-crosses/meta/tbl_sample_ped.txt
MINGQ=40

OUTDI=/data/anopheles/ag1000g/data/1000g_09_13/evaluation/phasing/phase1.AR2/callsets/autosome/truth/

# testing on pc
#PYBIN=/home/njh/PycharmProjects/phasing/bin
#OUTDI=/tmp/testing_pipe
#RAWH5=/home/njh/sync/3L_ag-cross.h5

WORKD=${OUTDI}/gq${MINGQ} 
FSTEM=3L_ag-cross.gq${MINGQ}

# This section makes the h5 file
mkdir $WORKD
$PYENV $PYBIN/filter_hdf5.py $RAWH5 $WORKD/$FSTEM -O -Q ${MINGQ} -P $PEDTB
touch ${WORKD}/${FSTEM}.h5.ok

$PYENV ${PYBIN}/hdf5_2_vcf.py ${WORKD}/${FSTEM}.h5 ${WORKD}/${FSTEM} -P $PEDTB
bgzip ${WORKD}/${FSTEM}.vcf
tabix -fp vcf ${WORKD}/${FSTEM}.vcf.gz
touch ${WORKD}/${FSTEM}.vcf.gz.ok
rm ${WORKD}/${FSTEM}.h5

# need to filter by the found mendelian errors
vcftools --gzvcf ${WORKD}/${FSTEM}.vcf.gz --exclude-positions \
  ${WORKD}/${FSTEM}_me.txt --recode --out ${WORKD}/${FSTEM}_me
mv ${WORKD}/${FSTEM}_me.recode.vcf ${WORKD}/${FSTEM}_me.vcf
bgzip ${WORKD}/${FSTEM}_me.vcf
tabix -fp vcf ${WORKD}/${FSTEM}_me.vcf.gz
touch ${WORKD}/${FSTEM}_me.vcf.gz.ok

# this section splits the vcf into the separate cross parts
$PYENV $PYBIN/split_vcf.py ${WORKD}/${FSTEM}_me.vcf.gz ${WORKD}/${FSTEM} \
  -P $PEDTB -B $VKEEP