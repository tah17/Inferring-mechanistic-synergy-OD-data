#!/bin/bash

jobID1=`qsub PBS/03g_all_drug_data_idxs.pbs`
jobID2=`qsub -W depend=afterok:$jobID1 PBS/03h_all_drug_fitting.pbs`
qsub -W depend=afterok:$jobID2 PBS/03i_all_drug_merge_fits.pbs