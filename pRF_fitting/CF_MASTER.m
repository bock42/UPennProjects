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

%% Make CF predictions
% outDir = fullfile('/data/jet/abock/cluster_shell_scripts/CF_scripts/AEK/10012014');
% session_dir = '/data/jet/abock/data/Template_Retinotopy/AEK/10012014/';
% subject_name = 'AEK_09242014_MPRAGE_ACPC_7T';
% runNums = [3,4,6];

% outDir = fullfile('/data/jet/abock/cluster_shell_scripts/CF_scripts/ASB/10272014');
% session_dir = '/data/jet/abock/data/Template_Retinotopy/ASB/10272014/';
% subject_name = 'ASB_10272014_MPRAGE_ACPC_7T';
% runNums = [2,4,6];

% outDir = fullfile('/data/jet/abock/cluster_shell_scripts/CF_scripts/GKA/10152014');
% session_dir = '/data/jet/abock/data/Template_Retinotopy/GKA/10152014/';
% subject_name = 'GKA_10152014_MPRAGE_ACPC_7T';
% runNums = [2,4,6];

% outDir = fullfile('/data/jet/abock/cluster_shell_scripts/CF_scripts/GKA/10152014_smooth_average');
% session_dir = '/data/jet/abock/data/Template_Retinotopy/GKA/10152014_smooth_average/';
% subject_name = 'GKA_10152014_MPRAGE_ACPC_7T';
% runNums = [2,4,6];

outDir = fullfile('/data/jet/abock/cluster_shell_scripts/CF_scripts/GKA/10152014_smooth_average/LGN/V1only');
session_dir = '/data/jet/abock/data/Template_Retinotopy/GKA/10152014_smooth_average/';
subject_name = 'GKA_10152014_MPRAGE_ACPC_7T';
runNums = [2,4,6];

% outDir = fullfile('/data/jet/abock/cluster_shell_scripts/CF_scripts/ASB/10272014_smooth_average/LGN/V1only');
% session_dir = '/data/jet/abock/data/Template_Retinotopy/ASB/10272014_smooth_average/';
% subject_name = 'ASB_10272014_MPRAGE_ACPC_7T';
% runNums = [2,4,6];

% outDir = fullfile('/data/jet/abock/cluster_shell_scripts/CF_scripts/AEK/10012014_smooth_average/LGN/V1only');
% session_dir = '/data/jet/abock/data/Template_Retinotopy/AEK/10012014_smooth_average/';
% subject_name = 'AEK_09242014_MPRAGE_ACPC_7T';
% runNums = [3,4,6];

templates = {'fine'};
hemis = {'lh' 'rh'};
srcROIs = {'LGN'};
srcfunc = 'dbrf.tf';
trgfunc = 's5.dbrf.tf';
cond = 'Movie';
V1only = 1;
cluster = 1;
create_CF_scripts(outDir,session_dir,subject_name,runNums,templates,hemis,...
    srcROIs,srcfunc,trgfunc,cond,V1only,cluster)
%% Plot and average maps
sessions = {...
    '/data/jet/abock/data/Template_Retinotopy/AEK/10012014_smooth_average/' ...
    '/data/jet/abock/data/Template_Retinotopy/ASB/10272014_smooth_average/' ...
    '/data/jet/abock/data/Template_Retinotopy/GKA/10152014_smooth_average/' ...
    };
subjects = {...
    'AEK_09242014_MPRAGE_ACPC_7T' ...
    'ASB_10272014_MPRAGE_ACPC_7T' ...
    'GKA_10152014_MPRAGE_ACPC_7T' ...
    };
subRuns = {[3,4,6],[2,4,6],[2,4,6]};
ROIs = {'LGN'};
hemis = {'lh' 'rh'};

template = 'fine';
srcfunc = 'dbrf.tf';
trgfunc = 's5.dbrf.tf';
cond = 'Movie';
DoG = 1;
V1only = 1; % V1 only!

for ss = 3%1:length(sessions)
    session_dir = sessions{ss};
    subject_name = subjects{ss};
    runNums = subRuns{ss};
    for roi = 1:length(ROIs)
        srcROI = ROIs{roi};
        for hh = 1:length(hemis)
            hemi = hemis{hh};
            for rr = 1:length(runNums)
                runNum = runNums(rr);
                plot_best_decimated_CF(session_dir,subject_name,runNum,hemi,...
                    srcROI,template,srcfunc,trgfunc,cond,DoG,V1only)
            end
        end
        average_CF(session_dir,subject_name,runNums,hemis,srcROI,template,srcfunc,trgfunc,cond,V1only)
    end
end
%% Create V1 corr mats
sessions = {...
    '/data/jet/abock/data/Template_Retinotopy/AEK/10012014_smooth_average' ...
    '/data/jet/abock/data/Template_Retinotopy/ASB/10272014_smooth_average' ...
    '/data/jet/abock/data/Template_Retinotopy/GKA/10152014_smooth_average' ...
    };
func = 's5.dbrf.tf';
areaFiles = {'lh.areas.vol' 'rh.areas.vol'};
hemis = {'lh' 'rh'};
regFile = 'brf_bbreg.dat';
allRuns = {[3,4,6] [2,4,6] [2,4,6]};
for i = 1:length(sessions)
    session_dir = sessions{i};
    b = find_bold(session_dir);
    runs = allRuns{i};
    for j = runs
        runDir = fullfile(session_dir,b{j});
        reg = fullfile(runDir,regFile);
        interp = 'nearest';
        inv = 1;
        for hh = 1:length(hemis)
            hemi = hemis{hh};
            outvol = fullfile(runDir,[areaFiles{hh} '.' func '.nii.gz']);
            invol = fullfile(runDir,[func '.nii.gz']);
            targvol = fullfile(session_dir,[areaFiles{hh} '.nii.gz']);
            mri_vol2vol(invol,targvol,outvol,reg,interp,inv)
            % Get average V1 signal
            funcAreas = load_nifti(fullfile(runDir,[areaFiles{hh} '.' func '.nii.gz']));
            V1ind = abs(funcAreas.vol) == 1;
            Funcvol = load_nifti(fullfile(runDir,[func '.nii.gz']));
            dims = size(Funcvol.vol);
            voltc = reshape(Funcvol.vol,dims(1)*dims(2)*dims(3),dims(4))';
            V1tc = mean(voltc(:,V1ind),2);
            V1corr = fisher_z_corr(corr(voltc,V1tc));
            % Save correlation matrix
            tmp = load_nifti(fullfile(runDir,[areaFiles{hh} '.' func '.nii.gz']));
            tmpcorr = reshape(V1corr,size(tmp.vol));
            tmp.vol = tmpcorr;
            save_nifti(tmp,fullfile(runDir,[hemi '.V1corr.' func '.nii.gz']));
            % Project to anat
            invol = fullfile(runDir,[hemi '.V1corr.' func '.nii.gz']);
            outvol = fullfile(runDir,[hemi '.V1corr.' func '.anat.nii.gz']);
            targvol = fullfile(session_dir,'lh.areas.vol.nii.gz');
            reg = fullfile(runDir,regFile);
            mri_vol2vol(invol,targvol,outvol,reg,interp);
        end
    end
end
%% Average runs
sessions = {...
    '/data/jet/abock/data/Template_Retinotopy/AEK/10012014_smooth_average' ...
    '/data/jet/abock/data/Template_Retinotopy/ASB/10272014_smooth_average' ...
    '/data/jet/abock/data/Template_Retinotopy/GKA/10152014_smooth_average' ...
    };
func = 's5.dbrf.tf';
areaFiles = {'lh.areas.vol' 'rh.areas.vol'};
hemis = {'lh' 'rh'};
regFile = 'brf_bbreg.dat';
allRuns = {[3,4,6] [2,4,6] [2,4,6]};
for i = 1:length(sessions)
    session_dir = sessions{i};
    for hh = 1:length(hemis)
        tmpvol = [];
        hemi = hemis{hh};
        ct = 0;
        b = find_bold(session_dir);
        runs = allRuns{i};
        for j = runs
            ct = ct + 1;
            runDir = fullfile(session_dir,b{j});
            tmp = load_nifti(fullfile(runDir,[hemi '.V1corr.' func '.anat.nii.gz']));
            tmpvol(ct,:,:,:) = tmp.vol;
        end
        avgV1corr = mean(tmpvol,1);
        save_nifti(tmp,fullfile(session_dir,[hemi '.V1corr.' func '.nii.gz']));
    end
end
%% Create LGN and MT masks
sessions = {...
    '/data/jet/abock/data/Template_Retinotopy/AEK/10012014_smooth_average/' ...
    '/data/jet/abock/data/Template_Retinotopy/ASB/10272014_smooth_average/' ...
    '/data/jet/abock/data/Template_Retinotopy/GKA/10152014_smooth_average/' ...
    };
subjects = {...
    'AEK_09242014_MPRAGE_ACPC_7T' ...
    'ASB_10272014_MPRAGE_ACPC_7T' ...
    'GKA_10152014_MPRAGE_ACPC_7T' ...
    };
hemis = {'lh' 'rh'};
maps = {'LGN' 'MT'};
SUBJECTS_DIR = getenv('SUBJECTS_DIR');
thresh = 5;
for ss = 1:length(sessions)
    session_dir = sessions{ss};
    subject_name = subjects{ss};
    for mm = 1:length(maps)
        map = maps{mm};
        for hh = 1:length(hemis)
            hemi = hemis{hh};
            in_vol = fullfile('~/data',[hemi '.' map '.cvs.nii.gz']);
            out_vol = fullfile(session_dir,[hemi '.' map '.prob.nii.gz']);
            ref_vol = fullfile(SUBJECTS_DIR,subject_name,'mri/T1.mgz');
            apply_cvs_inverse(in_vol,out_vol,ref_vol,subject_name,SUBJECTS_DIR);
            % Threshold voxels
            tmp = load_nifti(out_vol);
            tmp.vol(tmp.vol<thresh) = 0;
            tmp.vol(tmp.vol>0) = 1;
            save_nifti(tmp,fullfile(session_dir,[hemi '.' map '.nii.gz']));
            if strcmp(hemi,'rh')
                lh_vol_in = fullfile(session_dir,['lh.' map '.prob.nii.gz']);
                rh_vol_in = fullfile(session_dir,['rh.' map '.prob.nii.gz']);
                MNI_out = fullfile(session_dir,['mh.' map '.MNI.nii.gz']);
                mh_vol_out = fullfile(session_dir,['mh.' map '.prob.nii.gz']);
                map_type = 'foo'; % doesn't matter, as long as not 'copol' or 'varpol';
                average_xhemi_vol(session_dir,subject_name,lh_vol_in,rh_vol_in,MNI_out,mh_vol_out,map_type);
                tmp = load_nifti(mh_vol_out);
                tmp.vol(tmp.vol<thresh) = 0;
                tmp.vol(tmp.vol>0) = 1;
                save_nifti(tmp,fullfile(session_dir,['mh.' map '.nii.gz']));
            end
        end
    end
end
%% Average V1 correlation across runs
sessions = {...
    '/data/jet/abock/data/Template_Retinotopy/AEK/10012014_smooth_average/' ...
    '/data/jet/abock/data/Template_Retinotopy/ASB/10272014_smooth_average/' ...
    '/data/jet/abock/data/Template_Retinotopy/GKA/10152014_smooth_average/' ...
    };
subjects = {...
    'AEK_09242014_MPRAGE_ACPC_7T' ...
    'ASB_10272014_MPRAGE_ACPC_7T' ...
    'GKA_10152014_MPRAGE_ACPC_7T' ...
    };
func = 's5.dbrf.tf';
for ss = 1:length(sessions)
    session_dir = sessions{ss};
    subject_name = subjects{ss};
    lh_vol_in = fullfile(session_dir,['lh.V1corr.' func '.nii.gz']);
    rh_vol_in = fullfile(session_dir,['rh.V1corr.' func '.nii.gz']);
    MNI_out = fullfile(session_dir,['mh.V1corr.' func '.MNI.nii.gz']);
    mh_vol_out = fullfile(session_dir,['mh.V1corr.' func '.nii.gz']);
    map_type = 'foo'; % doesn't matter, as long as not 'copol' or 'varpol';
    average_xhemi_vol(session_dir,subject_name,lh_vol_in,rh_vol_in,MNI_out,mh_vol_out,map_type);
end
%% Define LGN ROI, using LGN mask and V1 correlation
sessions = {...
    '/data/jet/abock/data/Template_Retinotopy/AEK/10012014_smooth_average/' ...
    '/data/jet/abock/data/Template_Retinotopy/ASB/10272014_smooth_average/' ...
    '/data/jet/abock/data/Template_Retinotopy/GKA/10152014_smooth_average/' ...
    };
subjects = {...
    'AEK_09242014_MPRAGE_ACPC_7T' ...
    'ASB_10272014_MPRAGE_ACPC_7T' ...
    'GKA_10152014_MPRAGE_ACPC_7T' ...
    };
func = 's5.dbrf.tf';
thresh = 5;
LGNsize = 200;
for ss = 1:length(sessions)
    session_dir = sessions{ss};
    disp(session_dir);
    % load in average correlation
    V1corr = load_nifti(fullfile(session_dir,['mh.V1corr.' func '.nii.gz']));
    % Load LGN probabilistic mask
    for hh = 1:length(hemis)
        hemi = hemis{hh};
        LGN = load_nifti(fullfile(session_dir,[hemi '.LGN.prob.nii.gz']));
        LGNmask = LGN;
        LGNmask.vol = nan(size(LGN.vol));
        LGNind = find(LGN.vol>=thresh);
        % Find the top voxels using 'LGNsize'
        [corrVals,corrInd] = sort(V1corr.vol(LGNind));
        disp([hemi ' LGN thresh vals = ' num2str(corrVals(end-(LGNsize-1))) ' ' num2str(corrVals(end))]);
        LGNmask.vol(LGNind(corrInd(end-(LGNsize-1):end))) = 1;
        save_nifti(LGNmask,fullfile(session_dir,[hemi '.V1corr.' func '.LGN.nii.gz']));
    end
end
%% Save the LGN threshold values
% /data/jet/abock/data/Template_Retinotopy/AEK/10012014_smooth_average/
% lh LGN thresh vals = 0.19135 0.23689
% rh LGN thresh vals = 0.19173 0.25524
% /data/jet/abock/data/Template_Retinotopy/ASB/10272014_smooth_average/
% lh LGN thresh vals = 0.31237 0.46272
% rh LGN thresh vals = 0.28127 0.41453
% /data/jet/abock/data/Template_Retinotopy/GKA/10152014_smooth_average/
% lh LGN thresh vals = 0.2677 0.39629
% rh LGN thresh vals = 0.23842 0.37479
%% Project LGN to functional space
sessions = {...
    '/data/jet/abock/data/Template_Retinotopy/AEK/10012014_smooth_average/' ...
    '/data/jet/abock/data/Template_Retinotopy/ASB/10272014_smooth_average/' ...
    '/data/jet/abock/data/Template_Retinotopy/GKA/10152014_smooth_average/' ...
    };
hemis = {'lh' 'rh'};
func = 's5.dbrf.tf';
maps = {'LGN'};
regname = 'brf_bbreg.dat';
inv = 1; % inverse, so project anat to func
for ss = 1:length(sessions)
    session_dir = sessions{ss};
    d = find_bold(session_dir);
    for mm = 1:length(maps)
        map = maps{mm};
        for hh = 1:length(hemis)
            hemi = hemis{hh};
            for rr = 1:length(d);
                invol = fullfile(session_dir,[hemi '.V1corr.' func '.LGN.nii.gz']);
                targvol = fullfile(session_dir,d{rr},'single_TR.nii.gz');
                reg = fullfile(session_dir,d{rr},regname);
                outvol = fullfile(session_dir,d{rr},[hemi '.' map '.nii.gz']);
                mri_vol2vol(targvol,invol,outvol,reg,'nearest',inv);
            end
        end
    end
end
%%














%% View the output maps
% session_dir = '/data/jet/abock/data/Template_Retinotopy/AEK/10012014/';
% subject_name = 'AEK_09242014_MPRAGE_ACPC_7T';
% runNum = 2;

session_dir = '/data/jet/abock/data/Template_Retinotopy/ASB/10272014/';
subject_name = 'ASB_10272014_MPRAGE_ACPC_7T';
runNum = 4;

% session_dir = '/data/jet/abock/data/Template_Retinotopy/GKA/10152014/';
% subject_name = 'GKA_10152014_MPRAGE_ACPC_7T';
% runNum = 2;

hemi = 'lh';
srcROI = 'cortex';
template = 'pRF';

d = find_bold(session_dir);
R2 = load_nifti(fullfile(session_dir,'CFs',d{runNum},...
    [hemi '.' srcROI '.' template '.run' num2str(runNum) '.R2.cfs.nii.gz']));
ecc = load_nifti(fullfile(session_dir,'CFs',d{runNum},...
    [hemi '.' srcROI '.' template '.run' num2str(runNum) '.ecc.cfs.nii.gz']));
pol = load_nifti(fullfile(session_dir,'CFs',d{runNum},...
    [hemi '.' srcROI '.' template '.run' num2str(runNum) '.pol.cfs.nii.gz']));
areas = load_nifti(fullfile(session_dir,'pRFs',[template '_templates'],...
    [hemi '.areas.' template '.nii.gz']));
thresh = ~(abs(areas.vol)<=3) & R2.vol>0.15;
surface_plot('pol',pol.vol,subject_name,hemi,'inflated',thresh);
surface_plot('ecc',ecc.vol,subject_name,hemi,'inflated',thresh);

surface_plot('areas',areas.vol,subject_name,hemi,'inflated');

%%













%% Fix below
hemis = {'lh' 'rh'};
runs = [2 4 6];
decimation_level='0.1';
src_surf = 'inflated';
srcROI = 'volume';
for rr = runs
    for hh = 1:length(hemis)
        hemi = hemis{hh};
        plot_best_decimated_CF(session_dir,subject_name,rr,hemi,decimation_level,src_surf,srcROI);
        plot_best_V1_decimated_CF(session_dir,subject_name,rr,hemi,decimation_level,src_surf,srcROI);
    end
end
%% Average CF runs
%   Usage:
%   average_CF(session_dir,subject_name,runs,srcROI,trgROI,hemis)
map_type = 'movie';
runs = [2 4 6];
srcROI = 'volume';
template = 'fine';
average_CF(session_dir,subject_name,map_type,runs,srcROI,template);
average_V1_CF(session_dir,subject_name,map_type,runs,srcROI,template);


%% Deprecated








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

%%
compare_CF_pRF(session_dir,'cortex',template,map_type)


%% Project to cvs MNI space
% project_CF_cvsMNI(session_dir,subject_name);
%% Plot CF map
% plot_CF_maps(session_dir);




