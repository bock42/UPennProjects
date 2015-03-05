#!/bin/bash
$SGE_ROOT/bin/linux-x64/qsub -binding linear:12 -pe unihost 12 -l h_vmem=40.2G,s_vmem=40G -M mspits@sas.upenn.edu run_func_apply_reg_G092714A.sh
