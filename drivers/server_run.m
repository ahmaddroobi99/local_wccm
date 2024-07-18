clear

cases = loadcases();
casename = cases(3);

datarootfolder = 'data';

datafolder = fullfile(datarootfolder,casename);

pool = gcp('nocreate'); % if no pool, do not create new one.
if isempty(pool)
    pool = parpool();
end

%-------------------------------------------------------------------------%
% inputs
%-------------------------------------------------------------------------%


dataext = load(fullfile(datafolder, 'diagnostics.mat'));
fnamet = fullfile(datafolder, 't.bin'); fidt = fopen(fnamet,'r');
t = fread(fidt,'double'); fclose(fidt);
dt = t(2) - t(1);
nt = length(t);

pin.Tend = t(end);
pin.dt = dt;
pin.dtsave = dt*4;
pin.datafolder = datafolder;


qhc = coarsenqh(datafolder, pin.dtsave); % truth data


pin.fsavedata = true;

%-------------------------------------------------------------------------%
% damping
%-------------------------------------------------------------------------%

pin.damping_type = 'scaledep';  % d + mu |k|^2
pin.d = 0;
pin.mu = 0.05;

%-------------------------------------------------------------------------%
% spectrum
%-------------------------------------------------------------------------%

pin.E0 = 0.8;
pin.num_modes_flat = 0; % number of modes with flat spectrum

%-------------------------------------------------------------------------%
% analysis options
%-------------------------------------------------------------------------%

analysisopts.fsavefigs = true;
analysisopts.plot_times = floor([pin.Tend/4 pin.Tend/2]);


warning('off','MATLAB:prnRenderer:opengl') % turn off opengl nodisplay warning

clear F;

ntotal = 0;
for Nmodes = [8]
    for Ltracers = [10 20 40 80]
        for sigmax = [0.005, 0.01, 0.1, 0.2]
            
            obsdatafile = dir(fullfile(datafolder, sprintf('*Ltracers%i_sigmax%.2g.mat', Ltracers, sigmax)));
            if isempty(obsdatafile)
                error('obs data not found');
            end
            obsdata = load(fullfile(obsdatafile.folder,obsdatafile.name));
            
            
            pin.X = obsdata.X;
            pin.sigmax = obsdata.sigmax;
            pin.Ltracers = obsdata.Ltracers;
            
            
            fprintf('submitting job: ')
            fprintf('(2*N+1) = %i, Ltracers = %i, sigmax: %.2g\n', Nmodes, Ltracers, sigmax);
            
            
            pin.N = Nmodes;
            pin.spectrum_type = 'third';
            
            ntotal = ntotal + 1;
            F(ntotal) = parfeval(pool, @lagrangianfilter_baro, 3, pin);
            
        end
    end
end
fprintf('done submitting jobs...\n');



% fileid = fopen([pwd(),filesep,'runner.log'],'w');
% fprintf(fileid,'total number of cases %i\n', ntotal);
% result = cell(ntotal,1);
for ii = 1:ntotal
    [id, pin, m, Xt] = fetchNext(F);
%     [pin, m, Xt] = fetchOutputs(F(id));
    
%     try
%         analysis(pin, analysisopts, m, 0, qhc, Xt);
%     catch E
%         disp(E)
%     end

    duration = (F(id).FinishDateTime - F(id).StartDateTime);
    
    fprintf('finished case %i: (2*N+1) = %i, Ltracers = %i, sigmax: %.2g', id, pin.N, pin.Ltracers, pin.sigmax);
    fprintf('... total time %s... progress %.2f\n', string(duration), 100*ii/ntotal);
%     fprintf(fileid, strc);
end
% fclose(fileid);


%%





