% Load the data from data.mat
% loadedData = load('C:\Users\ahmad\OneDrive\Desktop\Local_baro\drivers\baroRK4_N128_beta0.05_d0.005_nu2.5e-28_forcing-sin_tau0-0.04_dt0.02_Nt100_IC-supplied_set-11-13-23-134958\N32_L10_sigmao0.1_Te2000_sp-third_da-unif_dt0.2\data.mat', 'pin', 'Xt', 'm');


%-------------------------------------------------------------------------%
% inputs
%-------------------------------------------------------------------------%
loadedData = load('N8_L10_sigmao0.1_Te2000_sp-kolmo_da-unif_dt0.2\hello_data.mat', 'm', 'Xt', 'pin');


% Extract the variables from the loaded data
pin = loadedData.pin; 
Xt = loadedData.Xt;
m = loadedData.m; 
Tend= pin.Tend;


damping_type = pin.damping_type;  % d + mu |k|^2
d = pin.d;
if ~isfield(pin,'mu')
    pin.mu = 0;
end
mu = pin.mu;

spectrum_type = pin.spectrum_type;
if ~isfield(pin,'E0')
    pin.E0 = 0;
end
if ~isfield(pin,'num_modes_flat')
    pin.num_modes_flat = 0;
end
E0 = pin.E0;
num_modes_flat = pin.num_modes_flat; % number of modes with flat spectrum


% sigmax = pin.sigmax;
sigmax = 0.1;
Ltracers = pin.Ltracers;
N = pin.N; % highest wavenumber
dt = pin.dt;
dtsave = pin.dtsave;

% defaults


if ~isfield(pin,'fsavedata')
    pin.fsavedata = false;
end
fsavedata = pin.fsavedata;


if ~isfield(pin,'datafolder')
    pin.datafolder = '';
end
datafolder = pin.datafolder;
%-------------------------------------------------------------------------%
% setup
%-------------------------------------------------------------------------%
F = 0; % zero forcing or periodic forcing only (not tested)

t = 0:dt:Tend+dt;
n_t = length(t);

dtsave = round(dtsave/dt)*dt; % round to a factor of dt
nsavef = dtsave/dt;
tsave = t(1:nsavef:length(t));
n_tsave = length(tsave);


fprintf('Total modes (2*N+1) = %i^2, max wavenumber resolved %d\n', N, N/2-1)
Ntotal = N*N;

L = 2*pi;
x = L/N*(0:N-1)';             % x coordinate
y = L/N*(0:N-1)';             % y coordinate
[xx,yy] = meshgrid(x,y);      % xy grid

kx = [0:N/2-1 0 -N/2+1:-1]; % x wave vector
ky = [0:N/2-1 0 -N/2+1:-1]; % y wave vector
[kxx, kyy] = meshgrid(kx, ky);
knorm = hypot(kxx,kyy);

projx = (1i*-kyy./knorm);
projy = (1i* kxx./knorm);
projx(zeromodes2d(N)) = 0;
projy(zeromodes2d(N)) = 0;

if strcmp(damping_type, 'scaledep')
    gamma = (d + mu*knorm.^2);
elseif strcmp(damping_type, 'unif')
    gamma = d*ones(N,N);
elseif strcmp(damping_type,'exp')
    gamma = d;
else
    error('unkown damping type')
end

if strcmp(spectrum_type, 'kolmo')
    sigma = sqrt(2*gamma.*(E0*knorm.^(-5/3))./knorm);
elseif strcmp(spectrum_type, 'com_kolmo')
    sigma = sqrt(2*gamma.*(E0*(knorm/num_modes_flat).^(-5/3)));
    sigma(knorm < num_modes_flat) = sqrt(2*gamma(knorm < num_modes_flat)*E0./knorm(knorm < num_modes_flat));
elseif strcmp(spectrum_type, 'third')
    sigma = sqrt(2*gamma.*(E0*knorm.^(-3))./knorm);
elseif strcmp(spectrum_type, 'com_third')
    sigma = sqrt(2*gamma.*(E0*(knorm/num_modes_flat).^(-3)));
    sigma(knorm < num_modes_flat) = sqrt(2*gamma(knorm < num_modes_flat)*E0./knorm(knorm < num_modes_flat));
elseif strcmp(spectrum_type, 'fourth')
    sigma = sqrt(2*gamma.*(E0*knorm.^(-4))./knorm);
elseif strcmp(spectrum_type, 'flat')
    sigma = sqrt(2*gamma.*(E0*knorm.^(0))./knorm);
elseif strcmp(spectrum_type, 'exp')
    if isscalar(pin.sigmax)
        sigma = pin.sigmax*ones(N,N);
    else
        sigma = pin.sigmax;
    end
else
    error('unkown spectrum type')
end

sigma(junkmodes2d(N)) = 0;
Sigma = sigma(:).*speye(Ntotal);
gamma(junkmodes2d(N)) = 0;
Gamma = gamma(:).' .* speye(Ntotal);


%-------------------------------------------------------------------------%
% save output parameter struct

pout.tsave = tsave;
pout.dtsave = dtsave;
pout.N = N;
pout.Ntotal = Ntotal;

pout.sigma = sigma;
pout.Sigma = Sigma;

pout.gamma = gamma;
pout.Gamma = Gamma;

pout.x = x;
pout.y = y;
pout.kx = kx;
pout.ky = ky;
pout.kxx = kxx;
pout.kyy = kyy;


pin.computed = pout; % save computed parameters

% Proceed with the existing variables and parameters
datafolder = pin.datafolder;

% Assuming the rest of the setup is correctly initialized
% Coarsen the QH data
qhc = coarsenqh(datafolder, pin.dtsave);
opts.renderer = 'zbuffer';

% Plotting options
opts.plot_times = [100, 200, 300];
opts.fsavefigs = true;

% Analyze the results
plotSingleTrialTimeSeries(Xt, m, 1, 2000);  % Example for the first trial and 2000th time point

% s_analysis(pin,m,Xt);
