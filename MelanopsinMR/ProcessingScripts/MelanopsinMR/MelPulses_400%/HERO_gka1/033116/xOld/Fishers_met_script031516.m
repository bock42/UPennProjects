% this script was used to convert F values to p values and then do Fisher's
% combined probability test for subject HERO_gka1, session 031516 under
% MaxMelPulse protocol.

%% Set up
session_dir = '/data/jag/MELA/HERO_gka1/031516';
subject_name = 'HERO_gka1_3T';
runNums = 1:9;
funcs = {...
   'FIR_raw'...
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