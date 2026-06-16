#!/bin/bash

jobID1=`qsub PBS/03a_no_drug_data_idxs.pbs`
jobID2=`qsub -W depend=afterok:$jobID1 PBS/03b_no_drug_fitting.pbs`
qsub -W depend=afterok:$jobID2 PBS/03c_no_drug_merge_fits.pbs