#!/bin/bash

jobID1=`qsub PBS/03d_indv_drug_data_idxs.pbs`
jobID2=`qsub -W depend=afterok:$jobID1 PBS/03e_indv_drug_fitting.pbs`
qsub -W depend=afterok:$jobID2 PBS/03f_indv_drug_merge_fits.pbs