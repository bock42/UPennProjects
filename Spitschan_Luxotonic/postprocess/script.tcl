sclv_read_binary_values 1
read_binary_curv; 
set curvflag 1; 
set forcegraycurvatureflag 1;
foreach file_name { cope1.feat_stats_zstat1.fsaverage_sym.lh.nii.gz cope1.feat_stats_zstat1.fsaverage_sym.rh.nii.gz cope1.feat_stats_zstat2.fsaverage_sym.lh.nii.gz cope1.feat_stats_zstat2.fsaverage_sym.rh.nii.gz cope1.feat_stats_zstat3.fsaverage_sym.lh.nii.gz cope1.feat_stats_zstat3.fsaverage_sym.rh.nii.gz cope2.feat_stats_zstat1.fsaverage_sym.lh.nii.gz cope2.feat_stats_zstat1.fsaverage_sym.rh.nii.gz cope2.feat_stats_zstat2.fsaverage_sym.lh.nii.gz cope2.feat_stats_zstat2.fsaverage_sym.rh.nii.gz cope2.feat_stats_zstat3.fsaverage_sym.lh.nii.gz cope2.feat_stats_zstat3.fsaverage_sym.rh.nii.gz cope3.feat_stats_zstat1.fsaverage_sym.lh.nii.gz cope3.feat_stats_zstat1.fsaverage_sym.rh.nii.gz cope3.feat_stats_zstat2.fsaverage_sym.lh.nii.gz cope3.feat_stats_zstat2.fsaverage_sym.rh.nii.gz cope3.feat_stats_zstat3.fsaverage_sym.lh.nii.gz cope3.feat_stats_zstat3.fsaverage_sym.rh.nii.gz } {

        set val $file_name
        sclv_read_from_dotw 0

        make_lateral_view
        redraw

        save_tiff ${file_name}-lateral.tiff


	make_lateral_view
	rotate_brain_y 180
	redraw

	   save_tiff ${file_name}-medial.tiff

}
exit
