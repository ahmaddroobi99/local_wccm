% clear;

cases = loadcases();
casename = cases(3);
display(casename);
datafolder = fullfile('data',casename);
display(datafolder);

dataext = load(fullfile(datafolder, 'diagnostics.mat'));
display(dataext);

fnamet = fullfile(datafolder, 't.bin'); fidt = fopen(fnamet,'r');
display(fnamet);

t = fread(fidt,'double'); fclose(fidt);
dt = t(2) - t(1);
nt = length(t);
Ntruth =  dataext.params.N;

%%%%%%%%%%%%%%%%%
% observations
%%%%%%%%%%%%%%%%%
Ltracers = 20;
sigmax = 0.20;
obsdatafile = dir(fullfile(datafolder, sprintf('obs_Ltracers%i_sigmax%.2g.mat', Ltracers, sigmax)));
if isempty(obsdatafile)
    error('obs data not found');
end
obsdata = load(fullfile(obsdatafile.folder,obsdatafile.name));


pin.Tend = t(end);

Nf = 16;
pin.Ltracers = Ltracers;
pin.N = Nf; % must be even, this is the total grid points in each dimension, e.g  16    32    64  
pin.X = obsdata.X;
pin.E0 = 0.8;
pin.damping_type = 'scale_dependent';
pin.d = 0;
pin.mu = 0.05;

pin.spectrum_type = 'third';
pin.num_modes_flat = 0;

pin.sigmax = obsdata.sigmax;
pin.Ltracers = obsdata.Ltracers;
pin.dt = dt;
pin.dtsave = dt*10;
pin.fsavedata = true;
pin.datafolder = casename;


paramstrain = train_msm(datafolder,dt*20);
params.train = 'msm';
d = (paramstrain.d - 1i*paramstrain.omega);
sigma = paramstrain.sigma;

k = [0:Ntruth/2 -Ntruth/2+1:-1]';
[KX, KY] = meshgrid(k,k);

k = [0:pin.N/2 -pin.N/2+1:-1]';
[KX, KY] = meshgrid(k,k);

pin.damping_type = 'scaledep';
pin.spectrum_type = 'kolmo';


[pin, m, Xt] = lagrangianfilter_baro(pin);

qhc = coarsenqh(datafolder,pin.dtsave);

opts.plot_times = [100 200];
opts.fsavefigs = true;
analysis(pin, opts, m, qhc, Xt)



