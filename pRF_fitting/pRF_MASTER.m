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
%% Project retinotopic templates to subject space
template_files = {...
    '~/data/2014-10-29.eccen-template.nii.gz' ...
    '~/data/2014-10-29.angle-template.nii.gz' ...
    '~/data/2014-10-29.areas-template.nii.gz'};
project_template(session_dir,subject_name,template_files,hemi)
%% Create occipital ROI
create_occipital(session_dir,subject_name);
%% xhemi check
% Checks that xhemireg and surfreg have been run for the specified
% freesurfer subject.
xhemi_check(session_dir,subject_name);
%% pRF analysis
% Generates a population receptive field (pRF) estimate using data obtained
%   while subjects viewed retinotopy stimuli (e.g. drifting bars). The
%   resulting pRF maps are then averaged across runs.
% If ROI = 'occipital', the averaged maps are plotted on the fsaverage_sym
%   surface, averaged across hemispheres, and converted to a format for
%   template fitting using Mathematica.
runs = [1,3,5];
ROI = 'cortex';
do_pRF(session_dir,subject_name,runs,ROI)
%% Run template fitting in Mathematica
% In a notebook in Mathematica, run the template fitting developed by Noah
%   Benson. The last line of that notebook saves a file (default:
%   ~/Desktop/template_fitting.mgz).
%% Create the pRF template .nii.gz files, using the output .mgz from Mathematica (above)
% Takes the .mgz output from Mathematica, convertes to nii.gz and separates
%   out the pol, ecc, and areas maps into individual volumes.
create_pRF_template(session_dir,subject_name)
%% Create the ccRF mat file
% Takes ~2 hours / run with a 12 core computer
% template = 'prf' % for cases
nruns = 6; % set the number of runs (in this case 6)
template = 'prf';
func = 'sdbrf.tf'
compute_distance = 1;
for rr = 1:nruns
    make_ccRF_mat(session_dir,subject_name,template,func,compute_distance)
end
%% ccRF analysis
% Takes ~2 hours / run with a 12 core computer
% Find cortico-cortical receptive fields (ccRF) within V1 for the specified
% roi <default - 'occipital'>
nruns = 6; % set the number of runs (in this case 6)
func = 'sdbrf.tf';
template = 'prf';
roi = 3; % 1 - V1; 2 - V1-V3 template; 3 - occipital; 4 - cortex; 5 - subcortical
for rr = 1:nruns
    do_ccRF(session_dir,template,rr,func,roi)
end
%% Plot ccRF maps based on roi (3 = occipital, 4 = cortex, 5 = subcortical)
nruns = 6; % set the number of runs (in this case 6)
template = 'prf';
func = 'sdbrf.tf';
roi = 3;
for rr = 1:nruns
    plot_ccRF(session_dir,subject_name,rr,func,template,roi)
end
%% Average ccRF maps across runs and hemispheres
template = 'prf';
func = 'sdbrf.tf';
condition = 'bars';
runs=[1 3 5]; % for AEK, bars = [1 2 5], movie = [3 4 6];
roi = 4;
average_ccRF_runs(session_dir,condition,runs,func,template,roi)