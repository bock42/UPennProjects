function outdir = func_apply_reg(session_dirs,func,out_dir)

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
outdir = fullfile(out_dir, 'timeseries');

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
        
        % Overwrite example_func2standard.nii.gz
        system(['flirt -in ' fullfile(session_dir,d{r}, func) ...
            ' -ref ' fullfile(session_dir,d{r},feat_dir,'reg','standard.nii.gz') ' -out ' ...
            fullfile(outdir, [func_new, '.timeseries.standard.nii.gz']) ...
            ' -init ' fullfile(session_dir,d{r},feat_dir,'reg','example_func2standard.mat') ...
            ' -applyxfm']);
        system(['flirt -in ' fullfile(session_dir,d{r}, 'firstlevel.feat', 'mean_func.nii.gz') ...
            ' -ref ' fullfile(session_dir,d{r},feat_dir,'reg','standard.nii.gz') ' -out ' ...
            fullfile(outdir, [func_new, '.mean.standard.nii.gz']) ...
            ' -init ' fullfile(session_dir,d{r},feat_dir,'reg','example_func2standard.mat') ...
            ' -applyxfm']);

    end
end
