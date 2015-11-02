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
%% Make contrasts
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
%% Create maps for subcortical ROIs
% AEK
%   LGNcenters = {[149 156 121] [107 157 123]}; % lh and rh
%   SCcenters = {[136 156 111] [123 152 114]}; % lh and rh
% ASB
%   LGNcenters = {[149 154 121] [106 153 121]}; % lh and rh
%   SCcenters = {[132 152 115] [123 152 115]}; % lh and rh
% GKA
%   LGNcenters = {[147 156 121] [107 155 124]}; % lh and rh
%   SCcenters = {[132 155 117] [123 153 118]}; % lh and rh
%
session_dir = '/jet/abock/data/Retinotopy/AEK/10012014';
feat_session_dir = '/jet/abock/data/Retinotopy/AEK/01152015';
LGNcenters = {[149 156 121] [107 157 123]}; % lh and rh
SCcenters = {[136 156 111] [123 152 114]}; % lh and rh

% session_dir = '/jet/abock/data/Retinotopy/ASB/10272014';
% feat_session_dir = '/jet/abock/data/Retinotopy/ASB/06022015';
% LGNcenters = {[149 154 121] [106 153 121]}; % lh and rh
% SCcenters = {[132 152 115] [123 152 115]}; % lh and rh
%
% session_dir = '/jet/abock/data/Retinotopy/GKA/06052015';
% feat_session_dir = '/jet/abock/data/Retinotopy/GKA/06052015';
% LGNcenters = {[147 156 121] [107 155 124]}; % lh and rh
% SCcenters = {[132 155 117] [123 153 118]}; % lh and rh

session_dir = '/jet/abock/data/Retinotopy/GKA/06052015';
feat_session_dir = '/jet/abock/data/Retinotopy/GKA/06052015';
LGNcenters = {[147 156 121] [107 155 124]}; % lh and rh
SCcenters = {[132 155 117] [123 153 118]}; % lh and rh

ROIs = {'LGN' 'SC'};
hemis = {'mh' 'lh' 'rh'};
srcROI = 'volume';
trgROI = 'prf_V1';
for rr = 1:length(ROIs);
    ROI = ROIs{rr};
    for hh = 1:length(hemis)
        hemi = hemis{hh};
        if strcmp(ROI,'LGN')
            if strcmp(hemi,'rh')
                center_vox = LGNcenters{2};
            else
                center_vox = LGNcenters{1};
            end
        elseif strcmp(ROI,'SC')
            if strcmp(hemi,'rh')
                center_vox = SCcenters{2};
            else
                center_vox = SCcenters{1};
            end
        end
        create_subcortical_CF_masks(session_dir,feat_session_dir,ROI,hemi,center_vox);
    end
end
%%











%% SC
input_vol = 'mh.volume.avg.co.prfs.nii.gz';
output_vol = 'lh.volume.SC.nii.gz';
center_vox = [123 154 120]; % GKA rh; [132 155 119]; % GKA lh; [122 151 116]; % ASB rh; [132 152 115]; % ASB lh; [121 152 117]; % AEK rh ; [134 153 116]; % AEK lh;
ROIsize = 75;
find_ROI_voxels(input_vol,output_vol,center_vox,ROIsize,voxsize);

%% LGN
input_vol = 'mh.volume.avg.co.prfs.nii.gz';
output_vol = 'lh.volume.LGN.nii.gz';
center_vox = [108 154 124]; % GKA rh; [146 155 120]; % GKA lh; [107 152 120]; % ASB rh; [149 152 118]; % ASB lh; [107 153 119]; % AEK rh ;[148 154 116]; % AEK lh ;
ROIsize = 200;
find_ROI_voxels(input_vol,output_vol,center_vox,ROIsize,voxsize);

%% Create subcortical ecc and pol pRF maps
ecc = load_nifti('./mh.volume.avg.coecc.prfs.nii.gz');
pol = load_nifti('./mh.volume.avg.copol.prfs.nii.gz');
lhLGN = load_nifti('./lh.volume.LGN.nii.gz');
lhSC = load_nifti('./lh.volume.SC.nii.gz');
rhLGN = load_nifti('./rh.volume.LGN.nii.gz');
rhSC = load_nifti('./rh.volume.SC.nii.gz');

nlhLGN = lhLGN;
nlhSC = lhSC;
nrhLGN = rhLGN;
nrhSC = rhSC;

nlhLGN.vol = ecc.vol .* lhLGN.vol;
save_nifti(nlhLGN,'./lh.volume.LGN.ecc.nii.gz');

nlhSC.vol = ecc.vol .* lhSC.vol;
save_nifti(nlhSC,'./lh.volume.SC.ecc.nii.gz');

nrhLGN.vol = ecc.vol .* rhLGN.vol;
save_nifti(nrhLGN,'./rh.volume.LGN.ecc.nii.gz');

nrhSC.vol = ecc.vol .* rhSC.vol;
save_nifti(nrhSC,'./rh.volume.SC.ecc.nii.gz');

nlhLGN.vol = pol.vol .* lhLGN.vol;
save_nifti(nlhLGN,'./lh.volume.LGN.pol.nii.gz');

nlhSC.vol = pol.vol .* lhSC.vol;
save_nifti(nlhSC,'./lh.volume.SC.pol.nii.gz');

nrhLGN.vol = pol.vol .* rhLGN.vol;
save_nifti(nrhLGN,'./rh.volume.LGN.pol.nii.gz');

nrhSC.vol = pol.vol .* rhSC.vol;
save_nifti(nrhSC,'./rh.volume.SC.pol.nii.gz');


%% Create subcortical ROIs
% Load in zstat volumes in fslview, create mask manually, thresholded at
% z > 2.58 (p < 0.01)
% Conservative estimate, to be included in SC and LGN, voxels needed a
% z-stat > 2.58 for BOTH the dots and checkerboard stimulus, for cope3/4
% (condition 1(2) > 2(1)).
%% xhemi check
% Check that xhemireg and
xhemi_check(session_dir,subject_name)
%% Project zstats to surface
session_dir = '/jet/abock/data/Retinotopy/ASB/06022015';
subject_name = 'ASB_10272014_MPRAGE_ACPC_7T';
run_dir = 'Series_007_BOLD_1.5mm_TR2000_mb4_P1_RUN9';
% subject_name = 'AEK_09242014_MPRAGE_ACPC_7T';
% session_dir = '/Users/abock/data/SC/AEK/01152015';
% run_dir = 'Series_008_bold_1.5_P2_mb4_TR2000_run3'; % run1 = dots, run2 = checker, run3 = M/P
feat_dir = 'sdbrf.tf.gfeat';
cope_dir = 'cope3.feat'; % condition 1 > condition 2
hemi = {'lh' 'rh'};
inname = fullfile(session_dir,run_dir,feat_dir,cope_dir,'stats','zstat1.nii.gz');
for h = 1:length(hemi);
    outname = fullfile(session_dir,[hemi{h} '_' run_dir '_' cope_dir '_zstat']);
    system(['mri_vol2surf --src ' inname ...
        ' --regheader ' subject_name ' --hemi ' hemi{h} ...
        ' --out ' outname '_surf.nii.gz --projfrac 0.5']);
    if h == 1;
        system(['mri_surf2surf --hemi ' hemi{h} ' --srcsubject ' ...
            subject_name ' --srcsurfval ' outname '_surf.nii.gz --trgsubject fsaverage_sym --trgsurfval ' ...
            outname '_fsavgsurf.nii.gz']);
    else
        system(['mri_surf2surf --hemi lh --srcsubject ' ...
            subject_name '/xhemi --srcsurfval ' outname '_surf.nii.gz --trgsubject fsaverage_sym --trgsurfval ' ...
            outname '_fsavgsurf.nii.gz']);
    end
    surface_plot('zstat',[outname '_fsavgsurf.nii.gz'],'fsaverage_sym');
end
%% zstat average
%session_dir = '/Users/abock/data/SC/AEK/01152015';
%run_dir = 'Series_008_bold_1.5_P2_mb4_TR2000_run3';
session_dir = '/Users/abock/data/SC/GKA/09242014';
run_dir = 'Series_007_BOLD_1.5mm_TR2000_mb4_P1_RUN9';
cope_dir = 'cope3.feat'; % 1 > 2
lhname = fullfile(session_dir,['lh_' run_dir '_' cope_dir '_zstat_fsavgsurf.nii.gz']);
rhname = fullfile(session_dir,['rh_' run_dir '_' cope_dir '_zstat_fsavgsurf.nii.gz']);
mhname = fullfile(session_dir,['mh_' run_dir '_' cope_dir '_zstat_fsavgsurf.nii.gz']);
lh = load_nifti(lhname);
rh = load_nifti(rhname);
mh = mean([lh.vol,rh.vol],2);
mh(abs(mh)<2.58 ) = nan; % 2.58 = p < 0.01; 1.96 = p < 0.05
lh.vol = mh;
save_nifti(lh,mhname);
surface_plot('zstat',mhname);
%% Hand draw label for MT and LOC on fsaverage_sym surface
% The above saves a 'mh_<run_dir>_<cope_dir>_zstat_fsavgsurf.nii.gz'
% overlay. Open this file in tksurfer on the fsaverage_sym surface, save as
% 'MT_fsavg_sym.label' and 'LOC_fsavg_sym.label' in session_dir
%% label 2 label
%session_dir = '/Users/abock/data/SC/AEK/01152015';
%subject_name = 'AEK_09242014_MPRAGE_ACPC_7T';
session_dir = '/Users/abock/data/SC/ASB/10012014';
subject_name = 'ASB_10012014_MPRAGE_ACPC_3T';
labels = {'MT' 'LOC'};
for l = 1:length(labels)
    srclabel = fullfile(session_dir,[labels{l} '_fsavg_sym.label']);
    lhtrglabel = fullfile(session_dir,['lh_' labels{l} '.label']);
    rhtrglabel = fullfile(session_dir,['rh_' labels{l} '.label']);
    system(['mri_label2label --srclabel ' srclabel ' --srcsubject fsaverage_sym' ...
        ' --trgsubject ' subject_name ' --trglabel ' lhtrglabel ...
        ' --regmethod surface --hemi lh']);
    system(['mri_label2label --srclabel ' srclabel ' --srcsubject fsaverage_sym' ...
        ' --trgsubject ' subject_name '/xhemi --trglabel ' rhtrglabel ...
        ' --regmethod surface --hemi lh']);
    copyfile(lhtrglabel,fullfile('/Applications/freesurfer/subjects/', ...
        subject_name,'label'));
    copyfile(rhtrglabel,fullfile('/Applications/freesurfer/subjects/', ...
        subject_name,'label'));
end
%% Project anatomical volumetric ROIs to functional space
% session_dir = '/Users/abock/data/Retinotopy/AEK/10012014';
% SC_dir = '/Users/abock/data/SC/AEK/01152015';
session_dir = '/Users/abock/data/Retinotopy/ASB/10272014';
SC_dir = '/Users/abock/data/SC/ASB/10012014';
hemi = {'lh' 'rh'};
subROIs = {'SC' 'LGN'};
func = 'sdbrf.tf';
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
for r = 1:nruns
    bbreg_out_file = fullfile(session_dir,d{r},'brf_bbreg.dat'); % registration file
    for ro = 1:length(subROIs)
        for h = 1:length(hemi)
            [~,~] = system(['mri_vol2vol --mov ' fullfile(session_dir,d{r},[func '.nii.gz']) ...
                ' --targ ' fullfile(SC_dir,[hemi{h} '_' subROIs{ro} '.nii.gz']) ' --o ' ...
                fullfile(session_dir,d{r},[hemi{h} '_' func '_' subROIs{ro} '.nii.gz']) ' --reg ' ...
                bbreg_out_file ' --inv --nearest']);
        end
    end
end
%% load in label vertices
%session_dir = '/Users/abock/data/Retinotopy/AEK/10012014';
%session_dir = '/Users/abock/data/SC/AEK/01152015';
%V1_session_dir = '/Users/abock/data/Retinotopy/AEK/10012014';
%subject_name = 'AEK_09242014_MPRAGE_ACPC_7T';
session_dir = '/Users/abock/data/Retinotopy/GKA/10152014';
V1_session_dir = '/Users/abock/data/Retinotopy/GKA/10152014';
subject_name = 'GKA_07302014_MPRAGE_ACPC_3T';
hemi = {'lh' 'rh'};
func = 'sdbrf.tf';
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
for h = 1:length(hemi)
    %%% Get surface ROI indices %%%
    V1 = load_nifti(fullfile(V1_session_dir,[hemi{h} '.areas_pRF.nii.gz']));
    ecc = load_nifti(fullfile(V1_session_dir,[hemi{h} '.ecc_pRF.nii.gz']));
    eval([hemi{h} '.V1 = find((V1.vol == -1 | V1.vol == 1) & ecc.vol <= 6.2116);']);% min visual angle radius for movie is 6.2116 (height)
    MT = read_label(subject_name,[hemi{h} '_MT']);
    eval([hemi{h} '.MT = MT(:,1) + 1;']);  % change from 0-based to 1-based index
    LOC = read_label(subject_name,[hemi{h} '_LOC']);
    eval([hemi{h} '.LOC = LOC(:,1) + 1;']); % change from 0-based to 1-based index
    for r = 1:nruns
        SC = load_nifti(fullfile(session_dir,d{r},[hemi{h} '_' func '_SC.nii.gz']));
        eval([hemi{h} '.SC_run' num2str(r) ' = find(SC.vol);']);
        LGN = load_nifti(fullfile(session_dir,d{r},[hemi{h} '_' func '_LGN.nii.gz']));
        eval([hemi{h} '.LGN_run' num2str(r) ' = find(LGN.vol);']);
    end
end
%% Compute average timecourse for each ROI, by run
%session_dir = '/Users/abock/data/Retinotopy/AEK/10012014';
%session_dir = '/Users/abock/data/SC/AEK/01152015';
session_dir = '/Users/abock/data/Retinotopy/ASB/10272014';
surfROIs = {'V1' 'MT' 'LOC'};
volROIs = {'SC' 'LGN'};
for h = 1:length(hemi)
    for r = 1:nruns
        surf = load_nifti(fullfile(session_dir,d{r},['sdbrf.tf_surf.' hemi{h} '.nii.gz']));
        surftc = squeeze(surf.vol);
        vol = load_nifti(fullfile(session_dir,d{r},'sdbrf.tf.nii.gz'));
        voltc = reshape(vol.vol,size(vol.vol,1)*size(vol.vol,2)*size(vol.vol,3),size(vol.vol,4));
        for su = 1:length(surfROIs)
            % Find ROI
            eval([hemi{h} '.' surfROIs{su} 'tc_run' num2str(r) ' = surftc(' hemi{h} '.' surfROIs{su} ',:);']);
            % Average across all vertices in ROI
            eval([hemi{h} '.' surfROIs{su} 'tc_run' num2str(r) ' = mean(' hemi{h} '.' surfROIs{su} 'tc_run' num2str(r) ',1);']);
        end
        for vo = 1:length(volROIs)
            % Find ROI
            eval([hemi{h} '.' volROIs{vo} 'tc_run' num2str(r) ' = voltc(' hemi{h} '.' volROIs{vo} '_run' num2str(r) ',:);']);
            % Average across all voxels in ROI
            eval([hemi{h} '.' volROIs{vo} 'tc_run' num2str(r) ' = mean(' hemi{h} '.' volROIs{vo} 'tc_run' num2str(r) ',1);']);
        end
    end
end
%% Convert to percent signal change
for h = 1:length(hemi)
    for r = 1:nruns
        for su = 1:length(surfROIs)
            eval(['tmp = ' hemi{h} '.' surfROIs{su} 'tc_run' num2str(r) ';']);
            dc = ones(size(tmp,1),1)*mean(tmp);
            tmp = ((tmp./dc) - 1) .*100;
            eval([hemi{h} '.' surfROIs{su} 'tc_run' num2str(r) ' = tmp;']);
        end
        for vo = 1:length(volROIs)
            eval(['tmp = ' hemi{h} '.' volROIs{vo} 'tc_run' num2str(r) ';']);
            dc = ones(size(tmp,1),1)*mean(tmp);
            tmp = ((tmp./dc) - 1) .*100;
            eval([hemi{h} '.' volROIs{vo} 'tc_run' num2str(r) ' = tmp;']);
        end
    end
end
%% Corrcoef - movie is runs 2,4,6 - ASB/GKA
clear X
X(:,1) = lh.LGNtc_run2;
X(:,2) = lh.SCtc_run2;
X(:,3) = lh.V1tc_run2;
X(:,4) = lh.MTtc_run2;
X(:,5) = lh.LOCtc_run2;
[lR1,lP1]=corrcoef(X);
clear X
X(:,1) = lh.LGNtc_run4;
X(:,2) = lh.SCtc_run4;
X(:,3) = lh.V1tc_run4;
X(:,4) = lh.MTtc_run4;
X(:,5) = lh.LOCtc_run4;
[lR2,lP2]=corrcoef(X);
clear X
X(:,1) = lh.LGNtc_run6;
X(:,2) = lh.SCtc_run6;
X(:,3) = lh.V1tc_run6;
X(:,4) = lh.MTtc_run6;
X(:,5) = lh.LOCtc_run6;
[lR3,lP3]=corrcoef(X);
clear X
X(:,1) = rh.LGNtc_run2;
X(:,2) = rh.SCtc_run2;
X(:,3) = rh.V1tc_run2;
X(:,4) = rh.MTtc_run2;
X(:,5) = rh.LOCtc_run2;
[rR1,rP1]=corrcoef(X);
clear X
X(:,1) = rh.LGNtc_run4;
X(:,2) = rh.SCtc_run4;
X(:,3) = rh.V1tc_run4;
X(:,4) = rh.MTtc_run4;
X(:,5) = rh.LOCtc_run4;
[rR2,rP2]=corrcoef(X);
clear X
X(:,1) = rh.LGNtc_run6;
X(:,2) = rh.SCtc_run6;
X(:,3) = rh.V1tc_run6;
X(:,4) = rh.MTtc_run6;
X(:,5) = rh.LOCtc_run6;
[rR3,rP3]=corrcoef(X);
%% Corrcoef - movie is runs 3,4,6 - AEK
clear X
X(:,1) = lh.LGNtc_run3;
X(:,2) = lh.SCtc_run3;
X(:,3) = lh.V1tc_run3;
X(:,4) = lh.MTtc_run3;
X(:,5) = lh.LOCtc_run3;
[R1,P1]=corrcoef(X);
clear X
X(:,1) = lh.LGNtc_run4;
X(:,2) = lh.SCtc_run4;
X(:,3) = lh.V1tc_run4;
X(:,4) = lh.MTtc_run4;
X(:,5) = lh.LOCtc_run4;
[R2,P2]=corrcoef(X);
clear X
X(:,1) = lh.LGNtc_run6;
X(:,2) = lh.SCtc_run6;
X(:,3) = lh.V1tc_run6;
X(:,4) = lh.MTtc_run6;
X(:,5) = lh.LOCtc_run6;
[R3,P3]=corrcoef(X);
%% Corrcoef - movie is runs 3,4,6
clear X
X(:,1) = rh.LGNtc_run3;
X(:,2) = rh.SCtc_run3;
X(:,3) = rh.V1tc_run3;
X(:,4) = rh.MTtc_run3;
X(:,5) = rh.LOCtc_run3;
[R1,P1]=corrcoef(X);
clear X
X(:,1) = rh.LGNtc_run4;
X(:,2) = rh.SCtc_run4;
X(:,3) = rh.V1tc_run4;
X(:,4) = rh.MTtc_run4;
X(:,5) = rh.LOCtc_run4;
[R2,P2]=corrcoef(X);
clear X
X(:,1) = rh.LGNtc_run6;
X(:,2) = rh.SCtc_run6;
X(:,3) = rh.V1tc_run6;
X(:,4) = rh.MTtc_run6;
X(:,5) = rh.LOCtc_run6;
[R3,P3]=corrcoef(X);
%% Corrcoef - dots are 1 4
clear X
X(:,1) = lh.LGNtc_run1;
X(:,2) = lh.SCtc_run1;
X(:,3) = lh.V1tc_run1;
X(:,4) = lh.MTtc_run1;
X(:,5) = lh.LOCtc_run1;
[R1,P1]=corrcoef(X);
clear X
X(:,1) = lh.LGNtc_run4;
X(:,2) = lh.SCtc_run4;
X(:,3) = lh.V1tc_run4;
X(:,4) = lh.MTtc_run4;
X(:,5) = lh.LOCtc_run4;
[R2,P2]=corrcoef(X);
%% Corrcoef - checker are 2 5
clear X
X(:,1) = lh.LGNtc_run2;
X(:,2) = lh.SCtc_run2;
X(:,3) = lh.V1tc_run2;
X(:,4) = lh.MTtc_run2;
X(:,5) = lh.LOCtc_run2;
[R1,P1]=corrcoef(X);
clear X
X(:,1) = lh.LGNtc_run5;
X(:,2) = lh.SCtc_run5;
X(:,3) = lh.V1tc_run5;
X(:,4) = lh.MTtc_run5;
X(:,5) = lh.LOCtc_run5;
[R2,P2]=corrcoef(X);
%% Corrcoef - MP are 3 6
clear X
X(:,1) = lh.LGNtc_run3;
X(:,2) = lh.SCtc_run3;
X(:,3) = lh.V1tc_run3;
X(:,4) = lh.MTtc_run3;
X(:,5) = lh.LOCtc_run3;
[R1,P1]=corrcoef(X);
clear X
X(:,1) = lh.LGNtc_run6;
X(:,2) = lh.SCtc_run6;
X(:,3) = lh.V1tc_run6;
X(:,4) = lh.MTtc_run6;
X(:,5) = lh.LOCtc_run6;
[R2,P2]=corrcoef(X);
%% Movie
lR = cat(3,lR1,lR2,lR3);
lP = cat(3,lP1,lP2,lP3);
lR = mean(lR,3);
lP = mean(lP,3);
lh.R = lR;
lh.zR = fisher_z_corr(lh.R);
rR = cat(3,rR1,rR2,rR3);
rP = cat(3,rP1,rP2,rP3);
rR = mean(rR,3);
rP = mean(rP,3);
rh.R = rR;
rh.zR = fisher_z_corr(rh.R);
%% Dots, Checker, M/P
R = cat(3,R1,R2);
P = cat(3,P1,P2);
R = mean(R,3);
P = mean(P,3);
%%
tmp = cat(3,lh.zR,rh.zR);
mh.zR = mean(tmp,3);

%%


figure;imagesc(lh.zR);colorbar
set(gca, 'XAxisLocation', 'top')
set(gca,'XTickLabel',{'','LGN','','SC','','V1','','MT','','LOC'},'FontSize',20)
set(gca,'YTickLabel',{'','LGN','','SC','','V1','','MT','','LOC'})
set(gca, 'Ticklength', [0 0])
caxis([0 1]);
colormap jet;
figure;imagesc(rh.zR);colorbar
set(gca, 'XAxisLocation', 'top')
set(gca,'XTickLabel',{'','LGN','','SC','','V1','','MT','','LOC'},'FontSize',20)
set(gca,'YTickLabel',{'','LGN','','SC','','V1','','MT','','LOC'})
set(gca, 'Ticklength', [0 0])
caxis([0 1]);
colormap jet;
figure;imagesc(mh.zR);colorbar
set(gca, 'XAxisLocation', 'top')
set(gca,'XTickLabel',{'','LGN','','SC','','V1','','MT','','LOC'},'FontSize',20)
set(gca,'YTickLabel',{'','LGN','','SC','','V1','','MT','','LOC'})
set(gca, 'Ticklength', [0 0])
caxis([0 1]);
colormap jet;

%%
figure;plot(lh.V1tc_run6);
figure;plot(lh.MTtc_run6);
figure;plot(lh.LOCtc_run6);
figure;plot(lh.SCtc_run6);
figure;plot(lh.LGNtc_run6);



















