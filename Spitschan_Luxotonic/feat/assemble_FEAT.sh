#!/bin/bash
subjID=$1

featBaseDir='/home/mspitschan/matlab/gkaguirrelab_Projects/Spitschan_Luxotonic/feat/'
dataBaseDir='/data/jag/mspitschan/Imaging/Protocols/MRLuxotonic'
cd $dataBaseDir/$subjID
dirs=`ls -d $PWD/RUN*`; 

mkdir $featBaseDir/$subjID
for f in $dirs; 
	do base=$(basename $f)
	echo $base
	sed 's,INPUT_DIR,'"$f"',g' $featBaseDir/design.fsf > $featBaseDir/$subjID/run_feat_$base.fsf
	echo feat $featBaseDir/$subjID/run_feat_$base.fsf > $featBaseDir/$subjID/run_feat_$base.sh
	echo $SGE_ROOT/bin/linux-x64/qsub -binding linear:8 -pe unihost 8 -l h_vmem=40.2G,s_vmem=40G -M mspits@sas.upenn.edu $featBaseDir/$subjID/run_feat_$base.sh >> $featBaseDir/$subjID/run_feat_submitjobs.sh 
done
