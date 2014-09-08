#! /bin/bash

/home/njh/pyenv/science/bin/python hdf5_2_tped.py \
	-i /data/anopheles/ag-crosses/data/release/0.1.GATK.PHASING.1000AG.AR2/h5/3L_ag-cross.h5 \
	-t /data/anopheles/ag1000g/data/1000g_09_13/evaluation/phasing/phase1.AR2/callsets/autosome/truth/samples.tfam \
	-r 3L:2000001-3000000 \
	-o /tmp/
