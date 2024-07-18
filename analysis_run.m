%-------------------------------------------------------------------------%
% load all .mat files in `datafolder` and run the analysis function on the
% datas
%-------------------------------------------------------------------------%

close all;
clear;



cases = loadcases();
casename = cases(3);
datarootfolder = 'data';
datafolder = fullfile(datarootfolder,casename);


dataext = load(fullfile(datafolder, 'diagnostics.mat'));
fnamet = fullfile(datafolder, 't.bin');
fidt = fopen(fnamet,'r'); t = fread(fidt,'double'); fclose(fidt);
dt = t(2) - t(1); nt = length(t);

files = dir(fullfile(datafolder, 'N*' ,'data.mat'));
analysisopts.fsavefigs = true;
analysisopts.plot_times = [t(end)/4,t(end)/2];


dtsave = dt*4; % INPUT
qhc = coarsenqh(datafolder, dtsave); % truth data

pool = gcp('nocreate'); % if no pool, do not create new one.
if isempty(pool)
    pool = parpool();
end
warning('off','MATLAB:prnRenderer:opengl') % turn off opengl nodisplay warning



T = table;
display(files);
% files= 'C:\Users\ahmad\OneDrive\Desktop\Research\LTM\barotropic-tracers-master\barotropic-tracers-master\data\baroRK4_N128_beta0.05_d0.005_nu2.5e-28_forcing-sin_tau0-0.04_dt0.02_Nt100_IC-supplied_set-11-13-23-134958\N16_L20_sigmao0.2_Te2_sp-kolmo_da-scaledep_dt0.2\data.mat'
parfor ii = 1:1
    
    ifile = files(ii);
    display(ifile);
    
    parts = strsplit(ifile.folder, filesep);
        display(parts);

    directory = parts{end};
            display(directory);

    fprintf('processing %s\n', directory)

    out = load(fullfile(ifile.folder, ifile.name));
    display(out.pin);
    disp( analysisopts);
    display(out.m);
    display(qhc);
    display(out.Xt);     
    Tnew = analysis(out.pin, analysisopts, out.m, qhc, out.Xt);
    T = [T; Tnew];

end


filterdiagfname = fullfile(datafolder,'filterdiagnostics.txt');
if exist(filterdiagfname,'file')
    delete(filterdiagfname);
end
writetable(T,filterdiagfname)


fprintf('done ... ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% run the script to add slopes on the Kolmogorov case for N=8 and plot
% streamfunction error

% figfolder = fullfile('figs',dataset); % INPUT
% 
% add_kolmogorovslopes
% 
% plot_streamfunctionerror


