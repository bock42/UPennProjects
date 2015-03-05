function outdir = func_apply_reg(session_dirs,func,out_dir)

%   Following feat_stats, this function will apply registrations so that
%   higher level feat will run.
%
%   inputs:
%   session_dirs - cell containing strings of the paths to the session
%   directories of interest
%   subject_name - freesurfer subject name
%   dirs - vector of feat directories (e.g. [1 4 7 10 13 16])
%   directories of interest
%   func - functional data file (default - 'dbrf.tf')
%
%   Written by Andrew S Bock Dec 2014
%
%   3/4/15  ms      Pulled out of feat_higher_level

%% Set default parameters
if ~exist('session_dirs','var')
    error('"session_dirs" not defined')
end
if ~exist('subject_name','var')
    error('"subject_name" not defined')
end
if ~exist('dirs','var')
    error('"dirs" not defined') % feat directories
end
if ~exist('func','var')
    func = 'dbrf.tf'; % functional data file
end
if ~exist('BACKUP', 'var')
    BACKUP = false;
end
if ~exist('OVERWRITE_REG_STANDARD', 'var')
    OVERWRITE_REG_STANDARD = false;
end

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
        % Overwrite example_func2standard.nii.gz
        system(['flirt -in ' fullfile(session_dir,d{r},feat_dir, func) ...
            ' -ref ' fullfile(session_dir,d{r},feat_dir,'reg','standard.nii.gz') ' -out ' ...
            fullfile(outdir, [func_new, '.standard.nii.gz') ...
            ' -init ' fullfile(session_dir,d{r},feat_dir,'reg','example_func2standard.mat') ...
            ' -applyxfm']);
    end
end