#!/bin/bash
$SGE_ROOT/bin/linux-x64/qsub -binding linear:2 -pe unihost 2 -l h_vmem=4.2G,s_vmem=4G -M mspits@sas.upenn.edu run_feat_first_level_apply_reg_A092714B.sh
