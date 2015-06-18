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
%% register_feat
% Registers the motion corrected and B0 unwarped functional volumes from
%   feat_mc_b0 to the corresponding Freesurfer anatomical image. See 'help
%   register_feat for details regarding default settings (e.g. despike,
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
%   the session_dir.
segment_anat(session_dir,subject_name);
%% Project anatomical ROIs to functional space
% Projects anatomical ROIs into functional space in each bold directory.
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
%   Usage:
%   run_CF(session_dir,subject_name,runNum,hemi,srcROI,trgROI,DoG,seedSig1,seedSig2,seedSig3,seedSig4)
run_CF(session_dir,subject_name);
%% Average CF runs
%   Usage:
%   average_CF(session_dir,subject_name,runs,srcROI,trgROI,hemis)
average_CF(session_dir,subject_name,runs);
%% Project to cvs MNI space
project_CF_cvsMNI(session_dir,subject_name);
%% Plot CF map
plot_CF_maps(session_dir);




