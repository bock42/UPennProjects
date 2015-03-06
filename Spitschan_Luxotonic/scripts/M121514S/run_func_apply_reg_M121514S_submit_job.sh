#!/bin/bash
$SGE_ROOT/bin/linux-x64/qsub -binding linear:6 -pe unihost 6 -l h_vmem=40.2G,s_vmem=40G -M mspits@sas.upenn.edu run_func_apply_reg_M121514S.sh
