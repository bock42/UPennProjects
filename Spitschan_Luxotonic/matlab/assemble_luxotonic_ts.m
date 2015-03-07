% tvvheFig = 1;
% nPlots = 3;
% plotInd = 1;
% roi = 'v1';
% session_dir = '/data/jag/mspitschan/Imaging/Protocols/MRLuxotonic/G092x14A/timeseries';
% luxotonic_ts_avg(session_dir, 'Isochromatic', theFig, nPlots, 1, roi)
%
% session_dir = '/data/jag/mspitschan/Imaging/Protocols/MRLuxotonic/A092x14B/timeseries';
% luxotonic_ts_avg(session_dir, 'Isochromatic', theFig, nPlots, 2, roi)
%
% session_dir = '/data/jag/mspitschan/Imaging/Protocols/MRLuxotonic/M12x514S/timeseries';
% luxotonic_ts_avg(session_dir, 'Isochromatic', theFig, nPlots, 3, roi)
close all;

theDirections = {'Isochromatic', 'LMS', 'Melanopsin'};

trSecs = 2;
nVols = 24;
t0 = (0:0.01:nVols-1)*trSecs;
t1 = (0:1:nVols-1)*trSecs;
close all
theFig = figure;
for d = 1:length(theDirections);
    
    dr = theDirections{d};
    
    nPlots = 4;
    plotInd = 1;
    roi = 'v1';
    session_dir = '/Data/jag/mspitschan/Imaging/Protocols/MRLuxotonic/G092x14A/timeseries';
    [t, cyc_avg1] = luxotonic_ts_avg(session_dir, dr, roi);
    subplot(4, 3, 4+(d-1));
    %plot(t0, 0.03*square(2*pi*1/48*t0)+0.2, '-k'); hold on;
    
    plot([min(t0) max(t0)], [0 0], '-k', 'Color', [0.5 0.5 0.5]); hold on;
    plot(t1, 100*mean(cyc_avg1, 2), '-k'); hold on;
    xlim([-1 48]); ylim([-0.31 0.31]);
    if d == 1
        ylabel({'Signal' 'change [%]'});
    end
    pbaspect([1 0.6 1]);
    set(0,'DefaultAxesTickDir', 'out'); box off; title('G092x14A');
    
    session_dir = '/Data/jag/mspitschan/Imaging/Protocols/MRLuxotonic/A092x14B/timeseries';
    [t, cyc_avg2] = luxotonic_ts_avg(session_dir, dr, roi);
    subplot(4, 3, 7+(d-1));
    %plot(t0, 0.03*square(2*pi*1/48*t0)+0.2, '-k'); hold on;
    
    plot([min(t0) max(t0)], [0 0], '-k', 'Color', [0.5 0.5 0.5]); hold on
    plot(t1, 100*mean(cyc_avg2, 2), '-k'); hold on;
    xlim([-1 48]); ylim([-0.31 0.31]);
    if d == 1
        ylabel({'Signal' 'change [%]'});
    end
    pbaspect([1 0.6 1]);
    set(0,'DefaultAxesTickDir', 'out'); box off;
    title('A092x14B');
    
    session_dir = '/Data/jag/mspitschan/Imaging/Protocols/MRLuxotonic/M12x514S/timeseries';
    [t, cyc_avg3] = luxotonic_ts_avg(session_dir, dr, roi);
    subplot(4, 3, 10+(d-1));
    %plot(t0, 0.03*square(2*pi*1/48*t0)+0.2, '-k'); hold on;
    plot([min(t0) max(t0)], [0 0], '-k', 'Color', [0.5 0.5 0.5]); hold on
    plot(t1, 100*mean(cyc_avg3, 2), '-k'); hold on;
    
    xlim([-1 48]); ylim([-0.31 0.31]); xlabel('Time [s]');
    if d == 1
        ylabel({'Signal' 'change [%]'});
    end
    pbaspect([1 0.6 1]);
    set(0,'DefaultAxesTickDir', 'out'); box off;
    title('M12x514S');w
    
    %group average
    subplot(4, 3, d);
    cyc_avg = [cyc_avg1  cyc_avg2 cyc_avg3];
    mean_cyc_avg = mean(cyc_avg, 2);
    sem = std(cyc_avg, [], 2)/sqrt(3);
    
    nVols =24;
    trSecs = 2;
    t0 = (0:0.01:nVols-1)*trSecs;
    plot([min(t0) max(t0)], [0 0], '-k', 'Color', [0.5 0.5 0.5]); hold on;
    plot(t0, 0.03*square(2*pi*1/48*t0+pi)-0.2, '-k'); hold on;
    shadedErrorBar(t1, 100*mean_cyc_avg, 100*sem, '-k');
    %plot([min(t0) max(t0)], [0 0], '-k', 'Color', [0.5 0.5 0.5]);
    
    xlim([-1 48]); ylim([-0.31 0.31]);
    if d == 1
        ylabel({'Signal' 'change [%]'});
    end
    pbaspect([1 0.6 1]);
    set(0,'DefaultAxesTickDir', 'out'); box off;
    title(dr);
    
    
end

subplot(4, 3, 1);     title(theDirections{1});
subplot(4, 3, 2);     title(theDirections{2});
subplot(4, 3, 3);     title(theDirections{3});
saveas(theFig, ['avg_' roi], 'pdf');
set(theFig, 'PaperPosition', [0 0 8 8]); %Position plot at left hand corner with width 12 and height 8.
set(theFig, 'PaperSize', [8 8]); %Set the paper to have width 12 and height 8.
