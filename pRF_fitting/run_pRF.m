session_dir = '/Users/abock/data/Retinotopy/ASB/10272014';
subject_name = 'ASB_10272014_MPRAGE_ACPC_7T';
runs = [1 3 5];
hemis = {'lh' 'rh'};
ROI = 'occipital';
d = listdir(fullfile(session_dir,'*bold_*'),'dirs');
out_dir = fullfile(session_dir,'prfs');
view_surfs = 1;
for r = 1:length(runs)
    run = runs(r);
    imFile = fullfile(session_dir,'Stimuli',['run' num2str(run)],'bars_images.mat');
    paramsFile = fullfile(session_dir,'Stimuli',['run' num2str(run)],'bars_params.mat');
    for h = 1:length(hemis)
        hemi = hemis{h};
        tcFile = fullfile(session_dir,d{r},['sdbrf.tf_surf.' hemi '.nii.gz']);
        pRF(session_dir,subject_name,run,hemi,tcFile,ROI,imFile,paramsFile,out_dir,view_surfs)
    end
end
        
        