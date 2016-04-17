%% This script will run the pRF analysis pipeline
%
%   *** It is assumed that the MRklar pipeline was already run ***
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

%% Copy template files
%   Get the template files from the Aguirre lab website
%       https://cfn.upenn.edu/aguirre/wiki/public:retinotopy_template
%
% copy to a local folder:
%
% template_files = {...
%   '~/data/2014-10-29.eccen-template.nii.gz' ...
%   '~/data/2014-10-29.angle-template.nii.gz' ...
%   '~/data/2014-10-29.areas-template.nii.gz'};
%% Project retinotopic templates to subject space
sessions = {...
    '/data/jet/abock/data/Template_Retinotopy/AEK/10012014/' ...
    '/data/jet/abock/data/Template_Retinotopy/ASB/10272014/' ...
    '/data/jet/abock/data/Template_Retinotopy/GKA/10152014/' ...
    };
subjects = {...
    'AEK_09242014_MPRAGE_ACPC_7T' ...
    'ASB_10272014_MPRAGE_ACPC_7T' ...
    'GKA_10152014_MPRAGE_ACPC_7T' ...
    };
for ss = 1:length(sessions)
    session_dir = sessions{ss};
    subject_name = subjects{ss};
    project_template(session_dir,subject_name)
end
%% Copy templates to 'pRFs/anat_templates
sessions = {...
    '/data/jet/abock/data/Template_Retinotopy/AEK/10012014/' ...
    '/data/jet/abock/data/Template_Retinotopy/ASB/10272014/' ...
    '/data/jet/abock/data/Template_Retinotopy/GKA/10152014/' ...
    };
subjects = {...
    'AEK_09242014_MPRAGE_ACPC_7T' ...
    'ASB_10272014_MPRAGE_ACPC_7T' ...
    'GKA_10152014_MPRAGE_ACPC_7T' ...
    };
hemis = {'lh' 'rh'};
templates = {'areas' 'ecc' 'pol'};
for ss = 1:length(sessions)
    session_dir = sessions{ss};
    subject_name = subjects{ss};
    outDir = fullfile(session_dir,'pRFs','anat_templates');
    if ~exist(outDir,'dir')
        mkdir(outDir);
    end
    for hh = 1:length(hemis)
        hemi = hemis{hh};
        for tt = 1:length(templates)
            template = templates{tt};
            system(['cp ' fullfile(session_dir,[hemi '.' template '.nii.gz']) ' ' ...
                fullfile(outDir,[hemi '.' template '.anat.nii.gz'])]);
        end
    end
end
%% Decimate the anatomical template files
sessions = {...
    '/data/jet/abock/data/Template_Retinotopy/AEK/10012014/' ...
    '/data/jet/abock/data/Template_Retinotopy/ASB/10272014/' ...
    '/data/jet/abock/data/Template_Retinotopy/GKA/10152014/' ...
    };
subjects = {...
    'AEK_09242014_MPRAGE_ACPC_7T' ...
    'ASB_10272014_MPRAGE_ACPC_7T' ...
    'GKA_10152014_MPRAGE_ACPC_7T' ...
    };
for ss = 1:length(sessions)
    session_dir = sessions{ss};
    subject_name = subjects{ss};
    tdir = fullfile(session_dir,'pRFs','anat_templates');
    decimate_templates(session_dir,subject_name,tdir);
end
%% pRF analysis
% Generates a population receptive field (pRF) estimate using data obtained
%   while subjects viewed retinotopy stimuli (e.g. drifting bars). The
%   resulting pRF maps are then averaged across runs.
% If ROI = 'occipital', the averaged maps are plotted on the fsaverage_sym
%   surface, averaged across hemispheres, and converted to a format for
%   template fitting using Mathematica.
run_pRF(session_dir,subject_name,runNum,hemi,srcROI)
%% Average pRF
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
allRuns = {[1,2,5],[1,3,5],[1,3,5]};
hemis = {'lh' 'rh'};
srcROI = 'LGN';
for ss = 1:length(sessions)
    session_dir = sessions{ss};
    subject_name = subjects{ss};
    runs = allRuns{ss};
    average_pRF(session_dir,subject_name,runs,srcROI);
end
%% Prepare for Mathematica
prepare_pRF_Mathematica(session_dir,subject_name)

%% Run template fitting in Mathematica
% In a notebook in Mathematica, run the template fitting developed by Noah
%   Benson.
%% Create the pRF template .nii.gz files, using the output .mgz from Mathematica (above)
% Takes the .mgz output from Mathematica, convertes to nii.gz and separates
%   out the pol, ecc, and areas maps into individual volumes.
sessions = {...
    '/data/jet/abock/data/Template_Retinotopy/AEK/10012014/' ...
    '/data/jet/abock/data/Template_Retinotopy/ASB/10272014/' ...
    '/data/jet/abock/data/Template_Retinotopy/GKA/10152014/' ...
    };
subjects = {...
    'AEK_09242014_MPRAGE_ACPC_7T' ...
    'ASB_10272014_MPRAGE_ACPC_7T' ...
    'GKA_10152014_MPRAGE_ACPC_7T' ...
    };
for ss = 1:length(sessions)
    session_dir = sessions{ss};
    subject_name = subjects{ss};
    create_pRF_template(session_dir,subject_name);
    tdir = fullfile(session_dir,'pRFs','pRF_templates');
    decimate_templates(session_dir,subject_name,tdir);
end
%% If doing correlation template fitting
%
%   see below
%
%%%
%% Create cluster shell scripts (pRF)
templateType = 'pRF';
sessions = {...
    '/data/jet/abock/data/Template_Retinotopy/AEK/10012014/' ...
    '/data/jet/abock/data/Template_Retinotopy/ASB/10272014/' ...
    '/data/jet/abock/data/Template_Retinotopy/GKA/10152014/' ...
    };
outDirs = {...
    ['/data/jet/abock/cluster_shell_scripts/fit_templates/AEK/10012014/' templateType] ...
    ['/data/jet/abock/cluster_shell_scripts/fit_templates/ASB/10272014/' templateType] ...
    ['/data/jet/abock/cluster_shell_scripts/fit_templates/GKA/10152014/' templateType] ...
    };
allRuns = {'[3,4,6]','[2,4,6]','[2,4,6]'};
func = 's5.dbrf.tf';
tcPart = 'full';
leaveOut = '0';
V2V3 = [0 1];
for ss = 1:length(sessions)
    session_dir = sessions{ss};
    for i = 1:length(V2V3)
        if V2V3(i)
            outDir = fullfile(outDirs{ss},'V2V3'); % V1-V2, V1-V3, AND V2-V3
            saveDir = fullfile(session_dir,'pRFs',templateType,func,'Movie','V2V3');
        else
            outDir = fullfile(outDirs{ss},'V1'); % V1-V2, V1-V3
            saveDir = fullfile(session_dir,'pRFs',templateType,func,'Movie','V1');
        end
        runs = allRuns{ss};
        create_regress_template_scripts(session_dir,templateType,outDir,runs,func,saveDir,tcPart,leaveOut,num2str(V2V3(i)));
    end
end
%% Create cluster shell scripts (anat)
templateType = 'anat';
sessions = {...
    '/data/jet/abock/data/Template_Retinotopy/AEK/10012014/' ...
    '/data/jet/abock/data/Template_Retinotopy/ASB/10272014/' ...
    '/data/jet/abock/data/Template_Retinotopy/GKA/10152014/' ...
    };
outDirs = {...
    ['/data/jet/abock/cluster_shell_scripts/fit_templates/AEK/10012014/' templateType] ...
    ['/data/jet/abock/cluster_shell_scripts/fit_templates/ASB/10272014/' templateType] ...
    ['/data/jet/abock/cluster_shell_scripts/fit_templates/GKA/10152014/' templateType] ...
    };
allRuns = {'[3,4,6]','[2,4,6]','[2,4,6]'};
func = 's5.dbrf.tf';
tcPart = 'full';
leaveOut = '0';
V2V3 = [0 1];
for ss = 1:length(sessions)
    session_dir = sessions{ss};
    for i = 1:length(V2V3)
        if V2V3(i)
            outDir = fullfile(outDirs{ss},'V2V3'); % V1-V2, V1-V3, AND V2-V3
            saveDir = fullfile(session_dir,'pRFs',templateType,func,'Movie','V2V3');
        else
            outDir = fullfile(outDirs{ss},'V1'); % V1-V2, V1-V3
            saveDir = fullfile(session_dir,'pRFs',templateType,func,'Movie','V1');
        end
        runs = allRuns{ss};
        create_regress_template_scripts(session_dir,templateType,outDir,runs,func,saveDir,tcPart,leaveOut,num2str(V2V3(i)));
    end
end
%% Convert coarse_model_templates 
%Takes the .mgz output from Mathematica, convertes to nii.gz and separates
%   out the pol, ecc, and areas maps into individual volumes.
%   note: this can also be run on the cluster
convert_Mathematica_templates(session_dir);

%% Decimate surfaces
% This has to be done in terminal IN LINUX!
%   cd $SUBJECTS_DIR/<subject_name>/surf
%   mris_decimate -d 0.1 ./lh.inflated ./lh.0.1.inflated
%   mris_decimate -d 0.1 ./rh.inflated ./rh.0.1.inflated

%% Decimate the pRF templates and bold runs
%   This can also be run on the cluster
tdir = fullfile(session_dir,'pRFs','coarse_model_templates');
decimate_templates(session_dir,subject_name,tdir);
decimate_bold(session_dir,subject_name,func);

%% Create cluster shell scripts (coarse)
templateType = 'coarse';
sessions = {...
    '/data/jet/abock/data/Template_Retinotopy/AEK/10012014/' ...
    '/data/jet/abock/data/Template_Retinotopy/ASB/10272014/' ...
    '/data/jet/abock/data/Template_Retinotopy/GKA/10152014/' ...
    };
outDirs = {...
    ['/data/jet/abock/cluster_shell_scripts/fit_templates/AEK/10012014/' templateType] ...
    ['/data/jet/abock/cluster_shell_scripts/fit_templates/ASB/10272014/' templateType] ...
    ['/data/jet/abock/cluster_shell_scripts/fit_templates/GKA/10152014/' templateType] ...
    };
allRuns = {'[3,4,6]','[2,4,6]','[2,4,6]'};
func = 's5.dbrf.tf';
tcPart = 'full';
leaveOut = '0';
V2V3 = [0 1];
for ss = 1:length(sessions)
    session_dir = sessions{ss};
    for i = 1:length(V2V3)
        if V2V3(i)
            outDir = fullfile(outDirs{ss},'V2V3'); % V1-V2, V1-V3, AND V2-V3
            saveDir = fullfile(session_dir,'pRFs',templateType,func,'Movie','V2V3');
        else
            outDir = fullfile(outDirs{ss},'V1'); % V1-V2, V1-V3
            saveDir = fullfile(session_dir,'pRFs',templateType,func,'Movie','V1');
        end
        runs = allRuns{ss};
        create_regress_template_scripts(session_dir,templateType,outDir,runs,func,saveDir,tcPart,leaveOut,num2str(V2V3(i)));
    end
end
%% Run 'coarse' template fits on the cluster ('regress_template')
% i.e. run the shell scripts created above

%% Find the best template
sessions = {...
    '/data/jet/abock/data/Template_Retinotopy/AEK/10012014/' ...
    '/data/jet/abock/data/Template_Retinotopy/ASB/10272014/' ...
    '/data/jet/abock/data/Template_Retinotopy/GKA/10152014/' ...
    };
hemis = {'lh' 'rh'};
templateType = 'coarse';
func = 's5.dbrf.tf';
fitType = 'V2V3'; % 'V1 = V1<->V2, V1<->V3; 'V2V3' =  V1<->V2, V1<->V3, AND V2<->V3
for ss = 1:length(sessions)
    session_dir = sessions{ss};
    disp(session_dir);
    for hh = 1:length(hemis)
        hemi = hemis{hh};
        disp(hemi);
        tdir = fullfile(session_dir,'pRFs',templateType,func,'Movie',fitType);
        [varexp,params,sorted_templates] = find_best_template(templateType,tdir,hemi,[],[],[],fitType);
        disp(params(1));
    end
end
%% Create fine templates, centered on the best template (above)
% In a notebook in Mathematica, create the fine templates

%% Convert fine_model_templates  
%Takes the .mgz output from Mathematica, convertes to nii.gz and separates
%   out the pol, ecc, and areas maps into individual volumes.
% note: this can also be run on the cluster
sessions = {...
    '/data/jet/abock/data/Template_Retinotopy/AEK/10012014/' ...
    '/data/jet/abock/data/Template_Retinotopy/ASB/10272014/' ...
    '/data/jet/abock/data/Template_Retinotopy/GKA/10152014/' ...
    };
subjects = {...
    'AEK_09242014_MPRAGE_ACPC_7T' ...
    'ASB_10272014_MPRAGE_ACPC_7T' ...
    'GKA_10152014_MPRAGE_ACPC_7T' ...
    };
V2V3 = [0 1];
for ss = 1:length(sessions)
    session_dir = sessions{ss};
    subject_name = subjects{ss};
    for i = 1:length(V2V3)
        if V2V3(i)
            tdir = fullfile(session_dir,'pRFs','fine_model_templates','V2V3');
        else
            tdir = fullfile(session_dir,'pRFs','fine_model_templates','V1');
        end
        convert_Mathematica_fine_templates(session_dir,tdir);
        decimate_templates(session_dir,subject_name,tdir);
    end
end
%% Create cluster shell scripts (fine)
templateType = 'fine';
sessions = {...
    '/data/jet/abock/data/Template_Retinotopy/AEK/10012014/' ...
    '/data/jet/abock/data/Template_Retinotopy/ASB/10272014/' ...
    '/data/jet/abock/data/Template_Retinotopy/GKA/10152014/' ...
    };
outDirs = {...
    ['/data/jet/abock/cluster_shell_scripts/fit_templates/AEK/10012014/' templateType] ...
    ['/data/jet/abock/cluster_shell_scripts/fit_templates/ASB/10272014/' templateType] ...
    ['/data/jet/abock/cluster_shell_scripts/fit_templates/GKA/10152014/' templateType] ...
    };
allRuns = {'[3,4,6]','[2,4,6]','[2,4,6]'};
func = 's5.dbrf.tf';
tcPart = 'full';
leaveOut = '0';
V2V3 = [0 1];
for ss = 1:length(sessions)
    session_dir = sessions{ss};
    for i = 1:length(V2V3)
        if V2V3(i)
            outDir = fullfile(outDirs{ss},'V2V3'); % V1-V2, V1-V3, AND V2-V3
            saveDir = fullfile(session_dir,'pRFs',templateType,func,'Movie','V2V3');
        else
            outDir = fullfile(outDirs{ss},'V1'); % V1-V2, V1-V3
            saveDir = fullfile(session_dir,'pRFs',templateType,func,'Movie','V1');
        end
        runs = allRuns{ss};
        create_regress_template_scripts(session_dir,templateType,outDir,runs,func,saveDir,tcPart,leaveOut,num2str(V2V3(i)));
    end
end
%% Run the 'fine' template fits on the cluster ('regress_template')
% i.e. run the shell scripts created above


