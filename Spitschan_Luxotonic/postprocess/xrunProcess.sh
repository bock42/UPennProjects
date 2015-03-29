cd $1/xrun.gfeat/surf

# First, run the tksurfer command that produces the images
tksurfer fsaverage_sym lh inflated  -tcl /home/mspitschan/matlab/gkaguirrelab_Projects/Spitschan_Luxotonic/postprocess/script.tcl

# Get rid of the margin
for f in `ls *.tiff`
do convert -shave 10%x22% $f $f
done

# Rename the TIFFs
rename -v -f  's/cope1/TransientON/g' *tiff
rename -v -f 's/cope2/TransientOFF/g' *tiff
rename -v -f 's/cope3/SustainedON/g' *tiff

rename -v -f 's/feat_stats_zstat1/Isochromatic/g' *tiff
rename -v -f 's/feat_stats_zstat2/LMS/g' *tiff
rename -v -f 's/feat_stats_zstat3/Melanopsin/g' *tiff

montage `ls TransientON.*` -tile 4x3 -geometry 250 TransientON.png
montage `ls TransientOFF.*` -tile 4x3 -geometry 250 TransientOFF.png
montage `ls SustainedON.*` -tile 4x3 -geometry 250 SustainedON.png

mkdir results
mv *png results
mv *tiff results


