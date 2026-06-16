#!/bin/bash

jobID1=`qsub PBS/s1a_prior_pred_data.pbs`
jobID2=`qsub -W depend=afterok:$jobID1 PBS/s1b_prior_pred_fit.pbs`
qsub -W depend=afterok:$jobID2 PBS/s1c_merge_fake_data_fits.pbs