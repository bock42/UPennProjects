%% Set up
session_dir = '/data/jag/MELA/HERO_aso1/033016';
subject_name = 'HERO_aso1_MaxMel';
runNums = 1:12;
funcs = {...
   'FIR_5mm' ...
    };
thresh = 0.05;
SUBJECTS_DIR = '/data/jag/MELA/freesurfer_subjects';

%% convert F values to p values
for ii = 1:length(funcs)
    func = funcs{ii}; 
convert_F_to_p(session_dir,subject_name,runNums,func,thresh,SUBJECTS_DIR)
end


%% Do the Fisher's combined probability test
for ii = 1:length(funcs)
    func = funcs{ii}; 
    fisher_combined_prob_test(session_dir,subject_name,runNums,func,thresh,SUBJECTS_DIR,true);
end