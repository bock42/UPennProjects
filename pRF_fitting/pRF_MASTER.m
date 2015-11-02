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
%% Create occipital ROI
create_occipital(session_dir,subject_name);
%% pRF analysis
% Generates a population receptive field (pRF) estimate using data obtained
%   while subjects viewed retinotopy stimuli (e.g. drifting bars). The
%   resulting pRF maps are then averaged across runs.
% If ROI = 'occipital', the averaged maps are plotted on the fsaverage_sym
%   surface, averaged across hemispheres, and converted to a format for
%   template fitting using Mathematica.
run_pRF(session_dir,subject_name,runNum,hemi,srcROI,imFileName,paramsFileName)
%% Average pRF
average_pRF(session_dir,subject_name,runs,srcROI);
%% Prepare for Mathematica
prepare_pRF_Mathematica(session_dir,subject_name)

%% Run template fitting in Mathematica
% In a notebook in Mathematica, run the template fitting developed by Noah
%   Benson.
%% Create the pRF template .nii.gz files, using the output .mgz from Mathematica (above)
% Takes the .mgz output from Mathematica, convertes to nii.gz and separates
%   out the pol, ecc, and areas maps into individual volumes.
create_pRF_template(session_dir,subject_name)
%% If doing correlation template fitting
%
%   see below
%
%%%
%% Takes the .mgz output from Mathematica, convertes to nii.gz and separates
%   out the pol, ecc, and areas maps into individual volumes.
% note: this can also be run on the cluster
convert_Mathematica_templates(session_dir);

%% Decimate surfaces
% This has to be done in terminal IN LINUX!
%   cd $SUBJECTS_DIR/<subject_name>/surf
%   mris_decimate -d 0.1 ./lh.inflated ./lh.0.1.inflated
%   mris_decimate -d 0.1 ./rh.inflated ./rh.0.1.inflated

%% Decimate the pRF templates and bold runs
%   This can also be run on the cluster
decimate_templates(session_dir,subject_name);
decimate_bold(session_dir,subject_name);

%% Create cluster shell scripts
%   session_dir = '/data/jet/abock/data/Retinotopy/AEK/10012014/';
%   outDir = '/data/jet/abock/cluster_shell_scripts/fit_templates/AEK';
%   runs = '[3,4,6]'; % must be a string (!)
create_regress_template_scripts(session_dir,outDir,runs);

%% Run template fits on the cluster ('regress_template')
% i.e. run the shell scripts created above
%% Find the best template
template = 'fine';
hemi = 'lh';
[varexp,params,sorted_templates] = find_best_template(session_dir,template,hemi);
%% Run template fine template fitting in Mathematica
% In a notebook in Mathematica, run the template fitting developed by Noah
%   Benson.
%% Takes the .mgz output from Mathematica, convertes to nii.gz and separates
%   out the pol, ecc, and areas maps into individual volumes.
% note: this can also be run on the cluster
convert_Mathematica_fine_templates(session_dir);
%% Decimate fine templates
decimate_fine_templates(session_dir,subject_name);
%% Create cluster shell scripts for fine template search
%   session_dir = '/data/jet/abock/data/Retinotopy/AEK/10012014/';
%   outDir = '/data/jet/abock/cluster_shell_scripts/fit_templates/AEK/fine_template_scripts';
%   runs = '[3,4,6]'; % must be a string (!)
create_regress_fine_template_scripts(session_dir,outDir,runs);

%% Plot the template fits as a 3D mesh
session_dirs = {...
    '/data/jet/abock/data/Retinotopy/AEK/10012014' ...
    '/data/jet/abock/data/Retinotopy/ASB/10272014' ...
    '/data/jet/abock/data/Retinotopy/GKA/10152014' ...
    };
template = 'fine';
plot_error_bars = 0;
plot_template_mesh(session_dirs,template,plot_error_bars);
%% Plot the best templates by template type
mat = plot_template_comparison;

%% Create cluster shell scripts for COARSE template search
%   session_dir = '/data/jet/abock/data/Retinotopy/AEK/10012014/';
%   outDir = '/data/jet/abock/cluster_shell_scripts/fit_templates/AEK/fine_template_scripts';
%   runs = '[3,4,6]'; % must be a string (!)
create_template_residual_scripts(session_dir,outDir,runs);

%% Create cluster shell scripts for FINE template search
%   session_dir = '/data/jet/abock/data/Retinotopy/AEK/10012014/';
%   outDir = '/data/jet/abock/cluster_shell_scripts/fit_templates/AEK/fine_template_scripts';
%   runs = '[3,4,6]'; % must be a string (!)
create_fine_template_residual_scripts(session_dir,outDir,runs);

%% Copy the best coase and fine template to the session_dir
templates  = {'coarse' 'fine'};
for tt = 1:length(templates)
    template = templates{tt};
    copy_best_template(session_dir,template);
end

%% Plot the template variance explained for each vertex
plot_template_varexp(session_dir,subject_name,hemi,template,makeplot)

%% Plot the template variance explained for each vertex
plot_fine_template_varexp(session_dir,subject_name,hemi,template,makeplot)
%% Deprecated



%% Plot the best template
% Makes a 3D scatter plot, and exports the variance explained,
%   corresponding template parameters, and distances from the best template
template = 'coarse';
hemi = 'rh';
makeplot=1;
[varexp,params,dists,newx,newy,newz] = plot_template_fits(session_dir,template,hemi,makeplot);
%% Make SC ROI
% Creates a sphere around the SC (left/right) in Freesurfer's cvs_MNI space,
%   then projects to subject native anatomical space
make_SC_sphere(session_dir,subject_name);
%% Make LGN ROI
% Takens an LGN ROI (left/right) in Freesurfer's cvs_MNI space, created
%   using FSL's 1mm atlas, and projects to subject native anatomical space
make_LGN_ROI(session_dir,subject_name);