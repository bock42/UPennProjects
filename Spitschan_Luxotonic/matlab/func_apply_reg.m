function out_dir = func_apply_reg(session_dirs,subjID,func,out_dir,calc_xrun)

%   Following feat_stats, this function will apply registrations on any
%   func file.
%
%   Written by Andrew S Bock Dec 2014
%
%   3/4/15  ms      Pulled out of feat_higher_level
%   3/5/15  ms      For use for dbrf.tf registration

%% Set default parameters
if ~exist('session_dirs','var')
    error('"session_dirs" not defined')
end
if ~exist('func','var')
    func = 'dbrf.tf'; % functional data file
end

feat_dir = 'firstlevel.feat';

mkdir(fullfile(out_dir, 'timeseries'));
out_dir_t = fullfile(out_dir, 'timeseries');

DO_VOL_REG = false;

%% Loop through session directories
for s = 1:length(session_dirs)
    session_dir = session_dirs{s};
    % Find bold run directories
    d = listdir(fullfile(session_dir,'*BOLD_*'),'dirs');
    if isempty(d)
        d = listdir(fullfile(session_dir,'*EPI_*'),'dirs');
    end
    if isempty(d)
        d = listdir(fullfile(session_dir,'RUN*'),'dirs');
    end
    nruns = length(d);
    % Copy over bbregister registration file
    % overwrite the example_func2standard.mat in each feat dir
    for r = 1:nruns
        
        % Figure out the name of the run. This is very specific to
        % Spitschan's naming convention.
        tmp = strsplit(d{r}, '_')
        
        direction = tmp{3};
        runNum = tmp{2};
        [~, id] = fileparts(session_dir);
        func_new = [direction '_' id '_' runNum];
        
        if DO_VOL_REG
            % First, convert per-run functional timeseries and mean to subject
            % anatomy
            system(['flirt -in ' fullfile(session_dir,d{r}, func) ...
                ' -ref ' fullfile(session_dir,d{r},feat_dir,'reg','standard.nii.gz') ' -out ' ...
                fullfile(out_dir, [func_new, '.timeseries.standard.nii.gz']) ...
                ' -init ' fullfile(session_dir,d{r},feat_dir,'reg','example_func2standard.mat') ...
                ' -applyxfm']);
            system(['flirt -in ' fullfile(session_dir,d{r}, 'firstlevel.feat', 'mean_func.nii.gz') ...
                ' -ref ' fullfile(session_dir,d{r},feat_dir,'reg','standard.nii.gz') ' -out ' ...
                fullfile(out_dir, [func_new, '.mean.standard.nii.gz']) ...
                ' -init ' fullfile(session_dir,d{r},feat_dir,'reg','example_func2standard.mat') ...
                ' -applyxfm']);
        end
        
        
%         
%         % Subtract mean and divide by mean
%         system(['fslmaths ' fullfile(session_dir,d{r}, func) ' -sub ' fullfile(session_dir,d{r}, 'firstlevel.feat', 'mean_func.nii.gz') ' -div ' fullfile(session_dir,d{r}, 'firstlevel.feat', 'mean_func.nii.gz') ' ' fullfile(out_dir_t, [func_new, '.timeseries.' subjID '_exf.nii.gz'])]);
%         
%         % Delete first 12 volumes
%         system(['fslroi ' fullfile(out_dir_t, [func_new, '.timeseries.' subjID '_exf.nii.gz']) ' ' fullfile(out_dir_t, [func_new, '.timeseries.' subjID '_exf.nii.gz']) ' 12 144']);
%         
%         
%         % Then, convert per-run functional time series and mean to the
%         % surface
%         % Time series
%         system(['mri_vol2surf --mov ' fullfile(out_dir_t, [func_new, '.timeseries.' subjID '_exf.nii.gz']) ' --reg ' fullfile(session_dir,d{r}, 'brf_bbreg.dat') ' --hemi lh --projfrac 0.5 --o ' fullfile(out_dir_t, [func_new, '.timeseries.' subjID '.lh.nii.gz'])]);
%         system(['mri_surf2surf --srcsubject ' subjID ' --sval ' fullfile(out_dir_t, [func_new, '.timeseries.' subjID '.lh.nii.gz']) ' --trgsubject fsaverage_sym --tval ' fullfile(out_dir_t, [func_new, '.timeseries.fsaverage_sym.lh.nii.gz']) ' --hemi lh']);
%         
        system(['mri_vol2surf --mov ' fullfile(out_dir_t, [func_new, '.timeseries.' subjID '_exf.nii.gz']) ' --reg ' fullfile(session_dir,d{r}, 'brf_bbreg.dat') ' --hemi rh --projfrac 0.5 --o ' fullfile(out_dir_t, [func_new, '.timeseries.' subjID '.rh.nii.gz'])]);
        system(['mri_surf2surf --srcsubject ' subjID '/xhemi --sval ' fullfile(out_dir_t, [func_new, '.timeseries.' subjID '.rh.nii.gz']) ' --trgsubject fsaverage_sym --tval ' fullfile(out_dir_t, [func_new, '.timeseries.fsaverage_sym.rh.nii.gz']) ' --hemi lh']);
        
    end
end

if calc_xrun
% Now, in the out dir, push some of the results through
theFunctionalFiles = {'cope1.feat/stats/zstat1.nii.gz' 'cope1.feat/stats/zstat2.nii.gz'  'cope1.feat/stats/zstat3.nii.gz'  ...
    'cope2.feat/stats/zstat1.nii.gz'  'cope2.feat/stats/zstat2.nii.gz'  'cope2.feat/stats/zstat3.nii.gz' ...
    'cope3.feat/stats/zstat1.nii.gz'  'cope3.feat/stats/zstat2.nii.gz'  'cope3.feat/stats/zstat3.nii.gz' };
mkdir(fullfile(out_dir, 'xrun.gfeat', 'surf'));
for f = 1:length(theFunctionalFiles)
    theFile = theFunctionalFiles{f};
    theFileNew = strrep(theFile, '/', '_');
    [~, tmp] = fileparts(theFileNew);
    [~, theFileNew] = fileparts(tmp);
    
%     system(['mri_vol2surf --mov ' fullfile(out_dir, 'xrun.gfeat', theFile) ' --regheader ' subjID ' --hemi lh --projfrac 0.5 --o ' fullfile(out_dir, 'xrun.gfeat', 'surf', [theFileNew, '.' subjID '.lh.nii.gz'])]);
%     system(['mri_surf2surf --srcsubject ' subjID ' --sval ' fullfile(out_dir, 'xrun.gfeat', 'surf', [theFileNew, '.' subjID '.lh.nii.gz']) ' --trgsubject fsaverage_sym --tval ' fullfile(out_dir, 'xrun.gfeat', 'surf', [theFileNew, '.fsaverage_sym.lh.nii.gz']) ' --hemi lh']);
%     
    system(['mri_vol2surf --mov ' fullfile(out_dir, 'xrun.gfeat', theFile) ' --regheader ' subjID ' --hemi rh --projfrac 0.5 --o ' fullfile(out_dir, 'xrun.gfeat', 'surf', [theFileNew, '.' subjID '.rh.nii.gz'])]);
    system(['mri_surf2surf --srcsubject ' subjID '/xhemi --sval ' fullfile(out_dir, 'xrun.gfeat', 'surf', [theFileNew, '.' subjID '.rh.nii.gz']) ' --trgsubject fsaverage_sym --tval ' fullfile(out_dir, 'xrun.gfeat', 'surf', [theFileNew, '.fsaverage_sym.rh.nii.gz']) ' --hemi lh']);
end
end
