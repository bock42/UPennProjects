%% This script will run a pRF analysis pipeline
%
% Required software:
%   Freesurfer, FSL
%
% Inputs:
%   session_dir - subject specific session directory
%       e.g. '~/data/Retinotopy/ASB/10012014/'
%   subject_name - freesurfer subject name
%       e.g. 'ASB_10012014_MPRAGE_ACPC_3T'
%
% Directory structure:
%   data_directory -> project directory -> subject directory ->
%       session directory -> dicom directory
%       e.g. ~/data/Retinotopy/ASB/10012014/DICOMS
%
% If physiological measures are collected, using the following directory
%   structure:
%   data_directory -> project directory -> subject directory ->
%       session directory -> physio directory
%       e.g. ~/data/Retinotopy/ASB/10012014/PulseOx
%
%   Written by Andrew S Bock Dec 2014

%% Sort_nifti
%   Sort dicoms into series specific directories, convert to nifti
sort_nifti(session_dir);
%% Freesurfer_7T
% If data were collected at 7T, see "Freesurfer_7T"
%% make_fieldmap
% Creates a B0 field map, and brain extracts the magnitude image by
%   registering that image to the freesurfer "brain.mgz" image.
%   If not already run, run subject through Freesurfer pipeline. "sort_nifti"
%   produced the typical commands for running through the Freesurfer pipeline
%   for data collect at 3T.
make_fieldmap(session_dir,subject_name);
%% skull_strip
% Creates skull stripped file MPRAGE_brain.nii.gz using FreeSurfer tools
skull_strip(session_dir,subject_name);
%% feat_mc_b0
% Motion correct and B0 unwarp functional runs. This script has the option
%   to despike data as well. The result will be a design file for feat, and
%   the corresponding script to run feat in terminal. See 'help feat_mc_b0'
%   for details regarding default settings (e.g. despike, TR, warp_dir).
feat_mc_b0(session_dir,subject_name);
%% bbregister
% Registers the motion corrected and B0 unwarped functional volumes from
%   feat_mc_b0 to the corresponding Freesurfer anatomical image. See 'help
%   bbregister for details regarding default settings (e.g. despike,
%   feat_dir, func).
register_feat(session_dir,subject_name);
%% Create regressors for denoise
% Creates a nuisance regressor text file, based on physiological noise
%   and motion. If a task-design, the motion parameters are made orthogonal
%   to the design.
create_regressors(session_dir);
%% Temporal Filter
% Remove temporal frequencies. Default is to use the 'detrend' function,
%   based on the 'type' input (see help temporal_filter). You can also pass
%   'bptf' as the 'type' input, which can run high, low, or band-pass
%   temporal filters.
temporal_filter(session_dir);
%% Segment freesurfer aseg.mgz volume
% Segments the freesurfer anatomical aseg.mgz volume into several ROIs in
%   the session_dir:
%
%   brain.nii.gz
%   aseg.gm.nii.gz
%   aseg.wm.nii.gz
%   aseg.lh_ventricle.nii.gz
%   aseg.rh_ventricle.nii.gz
%   aseg.third_ventricle.nii.gz
%   aseg.fourth_ventricle.nii.gz
%   aseg.brainstem.nii.gz
%   aseg.unknown.nii.gz
segment_anat(session_dir,subject_name);
%% Project anatomical ROIs to functional space
% Projects anatomical ROIs into functional space in each bold directory:
%   <func>.brain.nii.gz
%   <func>.aseg.gm.nii.gz
%   <func>.aseg.wm.nii.gz
%   <func>.aseg.lh_ventricle.nii.gz
%   <func>.aseg.rh_ventricle.nii.gz
%   <func>.aseg.third_ventricle.nii.gz
%   <func>.aseg.fourth_ventricle.nii.gz
%   <func>.aseg.brainstem.nii.gz
%   <func>.aseg.unknown.nii.gz
project_anat2func(session_dir,subject_name);
%% Create localWM timecourses
% Creates local white matter timecourses for each voxel, and saves these
%   timecourses as a 4D volume [func '.WMtc.nii.gz']. The default 'func'
%   input is 'brf', so the final 4D volume will be 'brf.WMtc.nii.gz'. Will
%   also return an output 4D matrix 'WMtc'.
create_localWMtc(session_dir);
%% Remove noise
% Removes physiological and other non-neuronal noise regressors
remove_noise(session_dir,subject_name);
%% Clean up
% Cleans up intermediate files and directories
clean_up(session_dir)
%% Smooth surface and volume
% Smooth the volume and/or surface functional volumes using a 5mm kernel
smooth_vol_surf(session_dir);
%% xhemi check
% Checks that xhemireg and surfreg have been run for the specified
% freesurfer subject.
xhemi_check(session_dir,subject_name);
%% Project retinotopic templates to subject space
project_template(session_dir,subject_name)
%% Create occipital ROI
create_occipital(session_dir,subject_name);

%% Calculate cortical distance (e.g. in V1)
hemis = {'lh' 'rh'};
for hh = 1:length(hemis)
    hemi = hemis{hh};
    calc_surface_distance(session_dir,subject_name,hemi);
end
%% run CF
run_CF(session_dir,subject_name)

%% Average CF runs
template = 'prf';
func = 'sdbrf.tf';
condition = 'movie';
runs=[2 4 6]; % for AEK, bars = [1 2 5], movie = [3 4 6];
roi = 5; % 1=V1; 2=V1_V3; 3=occipital; 4 = cortex; 5=subcortical;
hemi = {'lh'};
average_CF_runs(session_dir,condition,runs,func,template,roi,hemi)

%%
session_dir = '/Users/abock/data/Retinotopy/ASB/10272014';
subject_name = 'ASB_10272014_MPRAGE_ACPC_7T';
% seedSig1 = linspace(1,5,3); % millimeters
% seedSig2 = linspace(6,10,3);
% seedSig3 = linspace(0,10,3);
% seedSig = [seedSig1',seedSig2',seedSig3'];
% seedSig1 = (.5:.1:15)';
% seedSig2 = (1.25:0.25:2)';
% seedSig3 = (0:0.1:1)';

% Use this for actual analysis
% seedSig1 = (.5:.25:15)';
% seedSig2 = (1.25:0.25:3)';
% seedSig3 = (0:0.1:1)';
% Use this for testing
seedSig1 = (1:2:9)';
seedSig2 = (1:.5:3)';
seedSig3 = (0:0.25:1)';
seedSig4 = (0:0.25:1)';
seedSig = {seedSig1 seedSig2 seedSig3 seedSig4};

hemis = {'lh' 'rh'};
space = 'volume';
srcROI = 'volume';
trgROI = 'prf_V1';
DoG = 1; % difference of Gaussians
% Find bold run directories
d = listdir(fullfile(session_dir,'*BOLD_*'),'dirs');
if isempty(d)
    d = listdir(fullfile(session_dir,'*bold_*'),'dirs');
end
if isempty(d)
    d = listdir(fullfile(session_dir,'*EPI_*'),'dirs');
end
if isempty(d)
    d = listdir(fullfile(session_dir,'RUN*'),'dirs');
end
nruns = length(d);
disp(['Session_dir = ' session_dir]);
disp(['Number of runs = ' num2str(nruns)]);
for rr = [2 4 6];
    for hh = 1:length(hemis)
        hemi = hemis{hh};
        % Get source indices
        if strcmp(trgROI,'V1')
            areas = load_nifti(fullfile(session_dir,[hemi '.areas.nii.gz']));
        elseif strcmp(trgROI,'prf_V1')
            areas = load_nifti(fullfile(session_dir,[hemi '.areas_pRF.nii.gz']));
        end
        if strcmp(srcROI,'V3');
            V1ind = find(areas.vol<=1 & areas.vol >=-1);
            V1_3ind = find(areas.vol<=3 & areas.vol >=-3);
            V1_3ind(V1ind) = 0;
            srcind = find(V1_3ind);
        elseif strcmp(srcROI,'cortex')
            srcind = 1:length(areas.vol); % entire cortex
        elseif strcmp(srcROI,'volume')
            binfile = fullfile(session_dir,d{rr},'single_TR.nii.gz');
            src = load_nifti(binfile);
            srcind = find(src.vol > 0 & ~isnan(src.vol));
        end
        % Get target indices
        if strcmp(trgROI,'V1');
            V1 = load_nifti(fullfile(session_dir,[hemi '.areas.nii.gz']));
            trgind = find(V1.vol<=1 & V1.vol >=-1);
        elseif strcmp(trgROI,'prf_V1');
            V1 = load_nifti(fullfile(session_dir,[hemi '.areas_pRF.nii.gz']));
            trgind = find(V1.vol<=1 & V1.vol >=-1);
        end
        if strcmp(srcROI,'volume')
            srcfile = fullfile(session_dir,d{rr},'sdbrf.tf.nii.gz');
        elseif strcmp(srcROI,'cortex')
            srcfile = fullfile(session_dir,d{rr},['sdbrf.tf_surf.' hemi '.nii.gz']);
        end
        trgfile = fullfile(session_dir,d{rr},['sdbrf.tf_surf.' hemi '.nii.gz']);
        do_CF(session_dir,subject_name,rr,space,srcfile,trgfile,srcind,trgind,srcROI,trgROI,seedSig,hemi,DoG);
    end
end
%% Average CF maps across runs and hemispheres
template = 'prf';
func = 'sdbrf.tf';
condition = 'movie';
runs=[2 4 6]; % for AEK, bars = [1 2 5], movie = [3 4 6];
roi = 5; % 1=V1; 2=V1_V3; 3=occipital; 4 = cortex; 5=subcortical;
hemi = {'lh'};
average_CF_runs(session_dir,condition,runs,func,template,roi,hemi)