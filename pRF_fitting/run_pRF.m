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
%% pRF analysis
% Takes the denoised output from "denoise", projects the 4D timeseries to the surface, copies the Benson et al. (2013) template to subject surface space, creates an "occipital" ROI, and generates a population receptive field (pRF) estimate using stimuli generated from code based on the vista-disp directory of the mrVista toolbox. The resulting pRF maps from individual runs are averaged across runs, plotted on the fsaverage_sym surface, and then averaged across hemispheres. Finally, these maps are converted to a format for template fitting using Mathematica
make_pRF(session_dir,subject_name,runs)
        