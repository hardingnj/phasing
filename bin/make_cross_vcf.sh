#! /bin/bash
set -e
set -o pipefail

# This is a bash script that does much of the heavy lifting for the creation of VCF files
PYENV=/home/njh/pyenv/science/bin/python
PHBIN=/home/njh/git/phasing/bin
RAWH5=/data/anopheles/ag-crosses/data/release/0.1.GATK.PHASING.1000AG.AR2/h5/3L_ag-cross.h5 
PEDTB=/home/njh/git/ag-crosses/meta/tbl_sample_ped.txt
MINGQ=40
OUTDI=/data/anopheles/ag1000g/data/1000g_09_13/evaluation/phasing/phase1.AR2/callsets/autosome/truth/

# This section makes the h5 file
mkdir ${OUTDI}/gq${MINGQ} 
$PYENV $PYBIN/filter_hdf5.py $RAWH5 ${OUTDI}/gq${MINGQ}/3L_ag-cross.gq${MINGQ} -O -Q ${MINGQ} -P $PEDTB
touch ${OUTDI}/gq${MINGQ}/3L_ag-cross.gq${MINGQ}.h5.ok

$PYENV ${PYBIN}/hdf5_2_vcf.py ${OUTDI}/gq${MINGQ}/3L_ag-cross.gq${MINGQ}.h5 ${OUTDI}/gq${MINGQ}/3L_ag-cross.gq${MINGQ} -P $PEDTB
bgzip ${OUTDI}/gq${MINGQ}/3L_ag-cross.gq${MINGQ}.vcf
tabix -fp vcf ${OUTDI}/gq${MINGQ}/3L_ag-cross.gq${MINGQ}.vcf.gz
touch ${OUTDI}/gq${MINGQ}/3L_ag-cross.gq${MINGQ}.vcf.gz.ok
rm ${OUTDI}/gq${MINGQ}/3L_ag-cross.gq${MINGQ}.h5

# this section splits the vcf into the separate cross parts
# get samples 
vcfkeepsamples [SAMPLES] | perl -ne 'if ($_ =~ m/^#/ || $_ =~ m/0\/0|0\/1|1\/1/) { print $_;}' | bgzip -c >  ${CROSS}.vcf.gz


