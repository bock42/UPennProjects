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
%% Project pRF templates (if exist) to volume
session_dirs = {...
    '/data/jet/abock/data/Network_Connectivity/ASB/11042015' ...
    '/data/jet/abock/data/Network_Connectivity/GKA/11042015' ...
    };
subject_names = {...
    'A101415B' ...
    'G101415A' ...
    };
for ss = 1:length(session_dirs)
    session_dir = session_dirs{ss};
    subject_name = subject_names{ss};
    project_pRF2vol(session_dir,subject_name);
end
%% Make contrasts for dot stimulus
blockdur = 16.56;
make_covariates_dots_localizer(session_dir,runNum,blockdur)
%% for M vs P stimulus
blockdur = 12;
for runNum = [1 2 12 13]
    make_covariates_MP_localizer(session_dir,runNum,blockdur);
end
%% Create feat .fsf file to localize ROIs
% funcs = {...
%     'drf.tf' ...
%     's2.drf.tf' ...
%     's5.drf.tf' ...
%     };
% session_dirs = {...
%     '/data/jet/abock/data/Network_Connectivity/ASB/11042015' ...
%     '/data/jet/abock/data/Network_Connectivity/GKA/11042015' ...
%     };
% runNums = [1 2 12 13];
session_dirs = {...
    '/data/jet/abock/data/Retinotopy/ASB/10142015/' ...
    '/data/jet/abock/data/Retinotopy/GKA/10142015/' ...
    };
runNums = [3 4];
funcs = {...
    'rf.tf' ...
    'drf.tf' ...
    's2.rf.tf' ...
    's2.drf.tf' ...
    };
for ss = 1:length(session_dirs)
    session_dir = session_dirs{ss};
    for ff = 1:length(funcs)
        func = funcs{ff};
        for runNum = runNums
            feat_localizer(session_dir,runNum,func);
        end
    end
end
%% run feat
% run the feat scripts found in <session_dir>/feat_localizer_scripts

%% Run higher level (2 feat)
session_dirs = {...
    '/data/jet/abock/data/Retinotopy/ASB/10142015' ...
    '/data/jet/abock/data/Retinotopy/GKA/10142015' ...
    };
subject_names = {...
    'A101415B' ...
    'G101415A' ...
    };
funcs = {...
    's2.rf.tf' ...
    };
runNums = [3 4];
for ss = 1:length(session_dirs)
    session_dir = session_dirs{ss};
    subject_name = subject_names{ss};
    for ff = 1:length(funcs)
        func = funcs{ff};
        feat_higher_level_2feat_5copes(session_dir,subject_name,runNums,func);
    end
end
%% Run higher level (4 feat)
session_dirs = {...
    '/data/jet/abock/data/Network_Connectivity/ASB/11042015' ...
    '/data/jet/abock/data/Network_Connectivity/GKA/11042015' ...
    };
subject_names = {...
    'A101415B' ...
    'G101415A' ...
    };
funcs = {...
    'rf.tf' ...
    'drf.tf' ...
    's2.rf.tf' ...
    's2.drf.tf' ...
    };
runNums = [1 2 12 13];
for ss = 1:length(session_dirs)
    session_dir = session_dirs{ss};
    subject_name = subject_names{ss};
    for ff = 1:length(funcs)
        func = funcs{ff};
        feat_higher_level_4feat_5copes(session_dir,subject_name,runNums,func);
    end
end

%% Make LGN ROI
% Takens an LGN ROI (left/right) in Freesurfer's cvs_MNI space, created
%   using FSL's 1mm atlas, and projects to subject native anatomical space
session_dirs = {...
    '/data/jet/abock/data/Network_Connectivity/ASB/11042015' ...
    '/data/jet/abock/data/Network_Connectivity/GKA/11042015' ...
    };
session_dirs = {...
    '/data/jet/abock/data/Retinotopy/ASB/10142015' ...
    '/data/jet/abock/data/Retinotopy/GKA/10142015' ...
    };
subject_names = {...
    'A101415B' ...
    'G101415A' ...
    };
for ss = 1:length(session_dirs)
    session_dir = session_dirs{ss};
    subject_name = subject_names{ss};
    make_LGN_ROI(session_dir,subject_name);
end
%% Make SC ROI
% Creates a sphere around the SC (left/right) in Freesurfer's cvs_MNI space,
%   then projects to subject native anatomical space
make_SC_sphere(session_dir,subject_name);

%% Create LGN mask using anatomy and dot localizer
% runNum = higher level feat directory
session_dirs = {...
    '/data/jet/abock/data/Retinotopy/ASB/10142015' ...
    '/data/jet/abock/data/Retinotopy/GKA/10142015' ...
    };
funcs = {'s2.rf.tf' 's2.drf.tf'};
gfeatNum = 3;
for ss = 1:length(session_dirs)
    session_dir = session_dirs{ss};
    for ff = 1:length(funcs)
        func = funcs{ff};
        create_dot_LGN_mask(session_dir,func,gfeatNum);
    end
end
%% Copy above to any directory
%
%
%%
% session_dirs = {...
%     '/data/jet/abock/data/Network_Connectivity/ASB/11042015' ...
%     '/data/jet/abock/data/Network_Connectivity/GKA/11042015' ...
%     };
% gfeatNum = 1;
% funcs = {'rf.tf' 'drf.tf' 's2.rf.tf' 's2.drf.tf'};
% copes = 1:5;
% for ss = 1:length(session_dirs)
%     session_dir = session_dirs{ss};
%     for ff = 1:length(funcs)
%         func = funcs{ff};
%         convert_copes_to_psc(session_dir,func,gfeatNum,copes);
%     end
% end
%% Create M and P voxels with ROIs
session_dirs = {...
    '/data/jet/abock/data/Network_Connectivity/ASB/11042015' ...
    '/data/jet/abock/data/Network_Connectivity/GKA/11042015' ...
    };
subject_names = {...
    'A101415B' ...
    'G101415A' ...
    };
funcs = {'s2.rf.tf'};%{'rf.tf'};% 'drf.tf'  's2.drf.tf'};
gfeatNum = 1;
hemis = {'lh' 'rh'};
ROIs = {'LGN' 'V1_pRF' 'V2_pRF' 'V3_pRF'};
for ss = 1%:length(session_dirs)
    session_dir = session_dirs{ss};
    subject_name = subject_names{ss};
    for ff = 1:length(funcs)
        func = funcs{ff};
        for hh = 1:length(hemis)
            hemi = hemis{hh};
            for rr = 1:length(ROIs)
                ROI = ROIs{rr};
                create_M_vs_P_voxels(session_dir,subject_name,func,gfeatNum,hemi,ROI)
            end
        end
    end
end
%% Project ROIs to functional space
session_dirs = {...
    '/data/jet/abock/data/Network_Connectivity/ASB/11042015' ...
    '/data/jet/abock/data/Network_Connectivity/GKA/11042015' ...
    };
runNums = [3:11,14:22];
funcs = {'s2.rf.tf'};%{'rf.tf'};% 'drf.tf'  's2.drf.tf'};
for ss = 1%:length(session_dirs)
    session_dir = session_dirs{ss};
    for ff = 1:length(funcs)
        func = funcs{ff};
        project_MP_ROIs_2func(session_dir,runNums,func);
    end
end
%% Create cross correlation matrices
session_dirs = {...
    '/data/jet/abock/data/Network_Connectivity/ASB/11042015' ...
    '/data/jet/abock/data/Network_Connectivity/GKA/11042015' ...
    };
runNums = [3:11,14:22];
tcfuncs = {'s2.drf.tf'};
roifuncs = {'s2.rf.tf'};%{'rf.tf'};% 'drf.tf'  's2.drf.tf'};
for ss = 1%:length(session_dirs)
    session_dir = session_dirs{ss};
    for tf = 1:length(tcfuncs)
        tcfunc = tcfuncs{tf};
        for rf = 1:length(roifuncs)
            roifunc = roifuncs{rf};
            create_MP_correlation_matrix(session_dir,runNums,tcfunc,roifunc);
        end
    end
end
%%



%%% Deprecated %%%








%% Run an F-test
funcs = {...
    's2.drf.tf' ...
    's5.drf.tf' ...
    };
runs  = [3 4];
Comps = [1 2];
for ff = 1:length(funcs)
    func = funcs{ff};
    [F,df] = run_F_test(session_dir,subject_name,func,runs,Comps);
end
%% Make contrasts (old)
d = listdir(fullfile(session_dir,'*bold_*'),'dirs');
if isempty(d)
    d = listdir(fullfile(session_dir,'BOLD_*'),'dirs');
end
runNums = load(fullfile(session_dir,'runs.txt'));
blockdur = 16;
TR = 2;
for rr = 1:length(runNums);
    outputDir = fullfile(session_dir,d{rr});
    make_contrasts(runNums(rr),outputDir,blockdur,TR);
end
%% First level GLM (feat_stats)
func = 'sdbrf.tf';
% session_dirs = {...
%     '/Users/abock/data/SC/ASB/10012014/' ...
%     '/Users/abock/data/SC/ASB/10022014/' ...
%     '/Users/abock/data/SC/ASB/10082014/'};
%session_dirs = {'/jet/abock/data/Retinotopy/ASB/06022015'};
%session_dirs = {'/jet/abock/data/Retinotopy/GKA/06052015'};
session_dirs = {'/jet/abock/data/Retinotopy/AEK/01152015'};
for s = 1:length(session_dirs)
    session_dir = session_dirs{s};
    feat_stats(session_dir,func)
end
%% Higher level GLM (feat_higher_level)
% session_dirs = {...
%     '/Users/abock/data/SC/ASB/10012014/' ...
%     '/Users/abock/data/SC/ASB/10022014/' ...
%     '/Users/abock/data/SC/ASB/10082014/'};
%session_dirs = {'/jet/abock/data/Retinotopy/ASB/06022015'};
%subject_name = 'ASB_10272014_MPRAGE_ACPC_7T';
session_dirs = {'/jet/abock/data/Retinotopy/ASB/10012014'};
subject_name = 'ASB_10272014_MPRAGE_ACPC_7T';
% session_dirs = {'/jet/abock/data/Retinotopy/GKA/06052015'};
% subject_name = 'GKA_10152014_MPRAGE_ACPC_7T';
% session_dirs = {'/jet/abock/data/Retinotopy/AEK/01152015'};
% subject_name = 'AEK_09242014_MPRAGE_ACPC_7T';
%alldirs = {[1 4 7 10 13 16] [2 5 8 11 14 17] [3 6 9 12 15 18]};
%alldirs = {[1 3 5] [2 4 6]};
alldirs = {[1 2 3]};
%alldirs = {[1 4] [2 5] [3 6]};
%alldirs = {[1 3 5]};
%alldirs = {[1 2 4 5]};
func = 'sdbrf.tf';
design_file = '/Users/Shared/Matlab/bock42/MRI_preprocessing/feat_higher_level_template_3dirs.fsf';
%design_file = '/Users/Shared/Matlab/bock42/MRI_preprocessing/feat_higher_level_template_4dirs.fsf';
for dd = 1:length(alldirs)
    feat_higher_level(session_dirs,subject_name,alldirs{dd},func,design_file)
    %feat_higher_level(session_dirs,subject_name,alldirs{dd},func,design_file)
end
