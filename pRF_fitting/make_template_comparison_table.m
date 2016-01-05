function [sub_mat] = make_template_comparison_table

% Loads in subject data, creates comparison table

%% set defaults
subjects = {...
    'AEK_09242014_MPRAGE_ACPC_7T' ...
    'ASB_10272014_MPRAGE_ACPC_7T' ...
    'GKA_10152014_MPRAGE_ACPC_7T' ...
    };
session_dirs = {...
    '/data/jet/abock/data/Retinotopy/AEK/10012014' ...
    '/data/jet/abock/data/Retinotopy/ASB/10272014' ...
    '/data/jet/abock/data/Retinotopy/GKA/10152014' ...
    };
best_template_vals = {...
    [0.6 -0.2 -0.5; 0.0 -0.3 0.2] ...
    [0.4 -0.3 -0.3; 0.0 -0.2 0.2] ...
    [0.4 0.1 -0.4; 0.3 -0.1 -0.1] ...
    };
best_template_names = {...
    {'10.6.4' '4.5.11'} ...
    {'8.5.6' '4.6.11'} ...
    {'8.9.5' '7.7.8'} ...
    };
templates = {...
    'coarse' ...
    'fine' ...
    'pRF' ...
    'dpRF' ...
    'anat' ...
    };
hemis = {'lh' 'rh'};
%% Pull out data
sub_mat = zeros(length(session_dirs),length(hemis),length(templates));
progBar = ProgressBar(length(session_dirs),'Making matrix...');
for ss = 1:length(session_dirs)
    session_dir = session_dirs{ss};
    for hh = 1:length(hemis)
        hemi = hemis{hh};
        for tt = 1:length(templates)
            template = templates{tt};
            if strcmp(template,'coarse') || strcmp(template,'fine')
                [varexp] = ...
                    find_best_template(session_dir,template,hemi);
            elseif strcmp(template,'pRF')
                if strcmp(session_dir,'/data/jet/abock/data/Retinotopy/AEK/10012014')
                    runs = [3 4 6];
                else
                    runs = [2 4 6];
                end
                temp = best_template_names{ss}{hh};
                [tmpvarexp] = regress_template(session_dir,runs,hemi,temp,0);
                varexp = nansum(tmpvarexp(:,[2,3]),2);
            elseif strcmp(template,'dpRF')
                if strcmp(session_dir,'/data/jet/abock/data/Retinotopy/AEK/10012014')
                    runs = [3 4 6];
                else
                    runs = [2 4 6];
                end
                temp = 'pRF';
                [tmpvarexp] = regress_template(session_dir,runs,hemi,temp,0);
                varexp = nansum(tmpvarexp(:,[2,3]),2);
            elseif strcmp(template,'anat')
                if strcmp(session_dir,'/data/jet/abock/data/Retinotopy/AEK/10012014')
                    runs = [3 4 6];
                else
                    runs = [2 4 6];
                end
                temp = 'anat';
                [tmpvarexp] = regress_template(session_dir,runs,hemi,temp,0);
                varexp = nansum(tmpvarexp(:,[2,3]),2);
            end
            sub_mat(ss,hh,tt) = varexp(1);
        end
    end
    progBar(ss);
end