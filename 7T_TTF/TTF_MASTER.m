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
%clean_up(session_dir)
%% Smooth surface and volume
% Smooth the volume and/or surface functional volumes using a 5mm kernel
smooth_vol_surf(session_dir);
%% xhemi check
% Checks that xhemireg and surfreg have been run for the specified
% freesurfer subject.
xhemi_check(session_dir,subject_name);
%% Project retinotopic templates to subject space
project_template(session_dir,subject_name)

%% Make SC ROI
% Creates a sphere around the SC (left/right) in Freesurfer's cvs_MNI space,
%   then projects to subject native anatomical space
%   For GKA/08252015
%   center_voxs{1} = [132   155   119];
%   center_voxs{1} = [123   155   119];
%
make_SC_sphere(session_dir,subject_name,center_voxs)
%% Make LGN ROI
% Takens an LGN ROI (left/right) in Freesurfer's cvs_MNI space, created
%   using FSL's 1mm atlas, and projects to subject native anatomical space
make_LGN_ROI(session_dir,subject_name);

%% Make covariates for FSL's FEAT

make_covariates_7T_OneLight(session_dir,bold_num,mat_num,mat_dir)
%% Create feat script for running GLM
funcs = {...
    'dbrf.tf' ...
    'noWMdbrf.tf' ...
    };
for ff = 1:length(funcs)
    func = funcs{ff};
    feat_TTF(session_dir,func);
end
%% Run an F-test on bold runs from OneLight data at 7T
funcs = {...
    'rf.tf' ...
    };
for ff = 1:length(funcs)
    func = funcs{ff};
    [F,df] = run_F_test(session_dir,subject_name,func);
end
%% Combine hemispheres
combine_template_hemispheres(session_dir,subject_name,0);

%% Plot TTF
% creates (using make_TTF function) and plots TTFs for specified
% hemishperes and ROIs
% session_dir = '/data/jet/aguirre/Projects/Retinotopy/ASB/A102015B';
% subject_name = 'A101415B';
session_dir = '/data/jet/aguirre/Projects/Retinotopy/ASO/A102015O';
subject_name = 'A092115O';
hemi = 'mh';
template = 'anat';
func = 'rf.tf';
ROIs = {'V1low' 'V1mid' 'V1high'};%{'MT'};%{'V2/V3low' 'V2/V3mid' 'V2/V3high'}; % 'SC'; 'V1low'; 'V1mid'; 'V1high';
Fthresh = 4;
ctx_ROIs = [5 15 50]; % eccentricity bands
HzNames = {'Hz8','Hz10','Hz12','Hz16','Hz20','Hz24'};
% If entopic (otherwise use defaults)
xlims = [0 7];
ylims = [-1 5];
xTick = [0 1 2 3 4 5 6 7];
xLabels = {'0Hz','8Hz','10Hz','12Hz','16Hz','20Hz','24Hz','32Hz'};
for rr = 1:length(ROIs)
    ROI = ROIs{rr};
    [means{rr},sems{rr}] = psc_cope(session_dir,subject_name,hemi,template,func,ROI,Fthresh,ctx_ROIs,HzNames);
    plot_TTF(means{rr},sems{rr},ROI,xlims,ylims,xTick,xLabels);
end
%% Make main effect maps
% 'contrastNum' specifies the contrast for the main effect.
make_main_effect_map(session_dir,subject_name,func,contrastNum);



