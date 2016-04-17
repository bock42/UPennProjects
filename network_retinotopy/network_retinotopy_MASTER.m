%% This script will run a network retinotopy pipeline
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

%% Run standard preprocessing (e.g. see pRF_MASTER.m)
%   make note of the denoising steps used, we may wish to explore the
%   different choices at this stage.
%% Decimate the cortex so that V1 is ~50 vertices
%   This must be done in linux, as the 'mris_decimate' command does not
%   work on a Mac.
%
%   mris_decimate -d 0.01 ./lh.inflated ./lh.0.01.inflated

%% Decimate the timeseries
%   Obtain the timeseries data from each vertex in V1,V2,V3 (dorsal and
%   ventral)
session_dir = '/data/jet/abock/data/Template_Retinotopy/GKA/10152014/';
subject_name = 'GKA_10152014_MPRAGE_ACPC_7T';
func = 's5.dbrf.tf';
src_surf = 'inflated';
trg_surf = '0.01.inflated';
decimate_bold(session_dir,subject_name,func,src_surf,trg_surf);
%% Obtain the assigned visual areas, polar angles, eccentricities
%   Get these values for each vertex, using the pRF pipeline
subject_name = 'GKA_10152014_MPRAGE_ACPC_7T';
hemi = 'lh';
src_surf = 'inflated';
trg_surf = '0.01.inflated';
in_vol = '/data/jet/abock/data/Template_Retinotopy/GKA/10152014/pRFs/model_templates/lh.areas.7.9.6.nii.gz';
out_vol = '~/tmp.nii.gz';
decimate_surf(subject_name,hemi,src_surf,trg_surf,in_vol,out_vol);
surface_plot('blueareas',out_vol,subject_name,'lh',trg_surf);


