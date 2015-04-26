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
%% denoise
% Removes low frequencies, as well as non-neuronal signals, using pulseOx,
%   motion parameters, and anatomical ROIs (e.g. white matter, ventricles).
%   If a block design was used, ensures that motion parameter regressors
%   are orthogonal to the block contrasts.
denoise(session_dir,subject_name);
%% Clean up
% Cleans up intermediate files and directories
clean_up(session_dir)
%% Smooth surface and volume
% Smooth the volume and/or surface functional volumes using a 5mm kernel
func = 'dbrf.tf';
ROI = {'surface' 'volume'};
hemi = {'lh' 'rh'};
d = listdir(fullfile(session_dir,'*bold_*'),'dirs');
nruns = length(d);
poolobj = gcp; % Gets current pool, and if no pool, creates a new one
disp('Smoothing 4D timeseries...')
parfor rr = 1:nruns
    smooth_vol_surf(session_dir,rr,func,ROI,hemi)
end
delete(poolobj); % close parpool
disp('done.');
%% xhemi check
% Checks that xhemireg and surfreg have been run for the specified
% freesurfer subject.
xhemi_check(session_dir,subject_name);
%% Project retinotopic templates to subject space
template_files = {...
    '~/data/2014-10-29.eccen-template.nii.gz' ...
    '~/data/2014-10-29.angle-template.nii.gz' ...
    '~/data/2014-10-29.areas-template.nii.gz'};
project_template(session_dir,subject_name,template_files)
%% Create occipital ROI
create_occipital(session_dir,subject_name);

%% Calculate cortical distance (e.g. in V1)
ROI = 'prf_V1';
hemis = {'lh' 'rh'};
for hh = 1:length(hemis)
    hemi = hemis{hh};    
    if strcmp(ROI,'V1');
        V1 = load_nifti(fullfile(session_dir,[hemi '.areas.nii.gz']));
        ROIverts = find(V1.vol<=1 & V1.vol >=-1);
    elseif strcmp(ROI,'prf_V1');
        V1 = load_nifti(fullfile(session_dir,[hemi '.areas_pRF.nii.gz']));
        ROIverts = find(V1.vol<=1 & V1.vol >=-1);
    end
    calc_surface_distance(session_dir,subject_name,ROI,ROIverts,hemi);
end
%%
session_dir = '/Users/abock/data/Retinotopy/ASB/10272014';
subject_name = 'ASB_10272014_MPRAGE_ACPC_7T';
% seedSig1 = linspace(1,5,3); % millimeters
% seedSig2 = linspace(6,10,3);
% seedSig3 = linspace(0,10,3);
% seedSig = [seedSig1',seedSig2',seedSig3'];
seedSig1 = (.5:.5:15)';
seedSig2 = [1.1;2;4;8];
seedSig3 = (0:.5:2)';
seedSig = {seedSig1 seedSig2 seedSig3};
hemis = {'lh' 'rh'};
space = 'surface';
trgROI = 'V1'; %%%%%%% change this, after you re-calculate distance for prf_V1
func = 'sdbrf.tf';
DoG = 1; % difference of Gaussians
% Find bold run directories
d = listdir(fullfile(session_dir,'*BOLD_*'),'dirs');
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
        areas = load_nifti(fullfile(session_dir,[hemi '.areas.nii.gz'])); %%%%% adjust to pRF
        V1ind = areas.vol<=1 & areas.vol >=-1;
        V1_3ind = areas.vol<=3 & areas.vol >=-3;
        V1_3ind(V1ind) = 0;
        srcind = find(V1_3ind);
        %srcfile = fullfile(session_dir,d{rr},'sdbrf.tf.nii.gz');
        srcfile = fullfile(session_dir,d{rr},['sdbrf.tf_surf.' hemi '.nii.gz']);
        trgfile = fullfile(session_dir,d{rr},['sdbrf.tf_surf.' hemi '.nii.gz']);
        do_CF(session_dir,subject_name,rr,space,func,srcfile,trgfile,srcind,trgROI,seedSig,hemi,DoG);
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