
%FEATRegisterToFSAnat(0, 0, pwd,'G092014A-MRLuxotonic')
function [t, cyc_avg] = luxotonic_ts_avg(session_dir, direction, roi, minDeg, maxDeg)
%   Registers functional runs in feat to freesurfer anatomical
%
%   Usage:
%   register_functional_feat(session_dir,subject,func,SUBJECTS_DIR)
%
%   e.g. register_functional_feat(1,0,'~/data/ASB'),'ASB')
%
%   Defaults:
%   check = 0; do not check registration
%   adjust = 0; do not manually adjust registration
%   session_dir = NO DEFAULT, must define
%   subject = NO DEFAULT, must define
%   func = 'filtered_func_data'
%   hemi = {'lh','rh'};
%
%   Outputs:
%   The function will display the minimum cost value of registraion, for
%   each feat direcotry
%
%   The function will also create the following files:
%   # = number, for each stat file (e.g. cope1, tstat1, zstat1)
%   m = hemisphere, either l - left, or r - right
%   1) mh_surf_cope#.nii.gz; cope, unsmoothed
%   2) mh_smooth_surf_cope#.nii.gz; cope, smoothed on surface, 5mm kernel
%   3) mh_surf_tstat#.nii.gz; tstat, unsmoothed
%   4) mh_smooth_surf_tstat#.nii.gz; tstat, smoothed on surface, 5mm kernel
%   5) mh_surf_zstat#.nii.gz; zstatttt, unsmoothed
%   6) mh_smooth_surf_zstat#.nii.gz; zstat, smoothed on surface, 5mm kernel
%
%   Written by Andrew S Bock Sept 2014

%close all;
trSecs = 2;

%% Set up FSL variables
fsl_path = '/usr/local/fsl/';
setenv('FSLDIR',fsl_path)
setenv('FSLOUTPUTTYPE','NIFTI_GZ')
curpath = getenv('PATH');
setenv('PATH',sprintf('%s:%s',fullfile(fsl_path,'bin'),curpath));
SUBJECTS_DIR = getenv('SUBJECTS_DIR');
%% Find feat directories in the session directory

switch roi
    case 's1'
        s1_t = MRIread('/Volumes/PASSPORT/MRLuxotonic/Subjects/G092x14A/BOLD/label.s1.G092x14A.sym.lh.surf.mgh');
        s1 = find(s1_t.vol);
        x_roi = s1;
    case 'insula'
        a = dlmread('/Users/mspits/Desktop/insula.fsaverage_sym.label', ' ', 3);
        %s1_t = MRIread('/Volumes/PASSPORT/MRLuxotonic/Subjects/G092x14A/BOLD/label.insula.G092x14A.sym.lh.surf.mgh');
        %s1 = find(s1_t.vol);
        x_roi = a(:, 1);
        
    case 'v1'
        
        %% Get V1
        areas_t = MRIread(fullfile(SUBJECTS_DIR, 'fsaverage_sym/templates/2014-10-29.areas-template.nii.gz'));
        v1 = [find(areas_t.vol == 1) find(areas_t.vol == -1)];
        
        %% Get eccentricity map
        eccen_t = MRIread(fullfile(SUBJECTS_DIR, 'fsaverage_sym/templates/2014-10-29.eccen-template.nii.gz'));
        upperLimit = maxDeg;
        lowerLimit = minDeg;
        eccen_range = find(eccen_t.vol > lowerLimit & eccen_t.vol < upperLimit);
        
        %% Find the intersecting vertices
        x_roi = intersect(v1, eccen_range);
end

searchString = [direction '_*.brf.tf.timeseries*fsaverage_sym*'];

d = dir(fullfile(session_dir,searchString));
ts_v1_avg = [];



ts_v1_vox{length(x_roi)} = [];
for r = 1:length(d)
    % Get the time series
    ts = MRIread(fullfile(session_dir, d(r).name));
    ts = squeeze(ts.vol);
    ts_v1 = ts(x_roi, :);
    nVols = size(ts_v1, 2);
    ts_v1_avg = [ts_v1_avg mean(ts_v1)];
%     
%     % Experimental test
% for i = 1:size(ts_v1, 1);
%     ts_v1_vox{i} = [ts_v1_vox{i} mean(reshape(ts_v1(i, :), 24, 6), 2)];
% end
% 
%     
end
test_var = reshape(ts_v1_avg, 24, 144);

t1 = (0:1:nVols-1)*trSecs;
t = t1;
cyc_avg = mean(test_var, 2);



%
% % Make plots
% t0 = (0:0.01:nVols-1)*trSecs;
% t1 = (0:1:nVols-1)*trSecs;
% %title([direction '/V1 signal']);
% %subplot(2, 1, 1);
% % plot(t0, 0.03*square(2*pi*1/48*t0)+0.7, '-k'); hold on;
% % plot(t1, mean([ run(:).ts_v1_mean], 2)*100, '-k', 'LineWidth', 2);
% % plot([min(t0) max(t0)], [0 0], '--k');
% % xlim([0 312]); ylim([-0.9 0.9]); xlabel('Time [s]'); ylabel('BOLD Signal change [%]');
% % pbaspect([1 0.3 1]);
%
% %
% % % Look at cycle averages
% % xrun_ts_v1_mean = mean([ run(:).ts_v1_mean], 2);
% % xrun_ts_v1_mean(1:12) = [];
% % xrun_ts_v1_mean_per_cyc = reshape(xrun_ts_v1_mean, 24, 6);
% %
% % subplot(2, 1, 2);
% % shadedErrorBar([0:23]*2, 100*mean(xrun_ts_v1_mean_per_cyc, 2), 100*std(xrun_ts_v1_mean_per_cyc, [], 2)/sqrt(length(xrun_ts_v1_mean_per_cyc)))
% % xlim([0 47]); ylim([-0.4 0.4]); xlabel('Time [s]'); ylabel('BOLD Signal change [%]');
% % pbaspect([1 0.3 1]);
%
% tmp = mean([ run(:).ts_v1_mean], 2);
% % Get rid of first 12
% tmp(1:12) = [];
%
% %set(luxFig, 'PaperPosition', [0 0 20 12])
% %set(luxFig, 'PaperSize', [20 12]); %Set the paper to have width 5 and height 5.
% %saveas(luxFig, 'cycleAvg', 'pdf');
%
% figure(theFig)
% subplot(1, nPlots, plotInd);
% nCyc = 6;
% cyc_avg = 100*mean(reshape(tmp, 24, nCyc), 2);
% cyc_sem = 100*std(reshape(tmp, 24, nCyc), [], 2)/sqrt(nCyc);
%
% t = 0:2:47;
% t0 = 0:0.1:47;
% shadedErrorBar(t, cyc_avg, cyc_sem, {'LineWidth', 2, 'Color', 'k'}); hold on;
% plot(t0, 0.4+0.05*square(2*pi*1/48*t0-pi), '-k')
% plot([0 48], [0 0], '--k'); hold on;
% pbaspect([1 1 1]);
%
% switch roi
%     case 's1'
%             title([strrep(direction, '_', '')]);
%     case 'v1'
%         title([strrep(direction, '_', ''), ', ' num2str(lowerLimit) '-' num2str(upperLimit) ' deg']);
% end
%
% set(gca, 'XTick', [0 24 48]);
% set(gca, 'YTick', [-0.8 -0.6 -0.4 -0.2 0.0 0.2 0.4 0.6 0.8]);
% ylim([-0.55 0.55]);
% xlim([-1 48]);
%
%
% %% Save the figure
% set(gca,'TickDir','out')
% %set(gcf, 'PaperPosition', [0 0 12 4]);
% %set(gcf, 'PaperSize', [12 4]);
% %saveas(gcf, ['ts_v1_' direction], 'pdf');

