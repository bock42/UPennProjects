theFig = 1;
nPlots = 1;
plotInd = 1;
roi = 'v1';
session_dir = '/data/jag/mspitschan/Imaging/Protocols/MRLuxotonic/G092x14A/timeseries'; 
luxotonic_ts_avg(session_dir, 'Isochromatic', theFig, nPlots, 1, roi)

session_dir = '/data/jag/mspitschan/Imaging/Protocols/MRLuxotonic/A092x14B/timeseries'; 
luxotonic_ts_avg(session_dir, 'Isochromatic', theFig, nPlots, 2, roi)

session_dir = '/data/jag/mspitschan/Imaging/Protocols/MRLuxotonic/M12x514S/timeseries'; 
luxotonic_ts_avg(session_dir, 'Isochromatic', theFig, nPlots, 3, roi)