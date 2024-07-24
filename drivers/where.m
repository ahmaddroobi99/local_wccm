function [pin, m, Xt] = lagrangianfilter_baro_EnKF(pin)

%-------------------------------------------------------------------------%
% inputs
%-------------------------------------------------------------------------%

Tend = pin.Tend;
damping_type = pin.damping_type;  % d + mu |k|^2
d = pin.d;
if ~isfield(pin,'mu')
    pin.mu = 0;
end
mu = pin.mu;
spectrum_type = pin.spectrum_type;
E0 = pin.E0;
num_modes_flat = pin.num_modes_flat; % number of modes with flat spectrum
sigmax = pin.sigmax;
Ltracers = pin.Ltracers;
N = pin.N; % highest wavenumber
dt = pin.dt;
dtsave = pin.dtsave;
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
    if isscalar(pin.sigma)
        sigma = pin.sigma*ones(N,N);
    else
        sigma = pin.sigma;
    end
else
    error('unkown spectrum type')
end

sigma(junkmodes2d(N)) = 0;
Sigma = sigma(:).*speye(Ntotal);
gamma(junkmodes2d(N)) = 0;
Gamma = gamma(:).*speye(Ntotal);

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

%-------------------------------------------------------------------------%
% solve filter for velocity recovery (EnKF)
%-------------------------------------------------------------------------%

% R = zeros(N*N, N*N, n_tsave);
m = zeros(N*N, n_tsave);
X = cell2mat(pin.X);
dXdt = [(X(:,2)-X(:,1)) (X(:,2:end) - X(:,1:end-1))]./dt;

% Initialize ensembles
num_ensembles = 15;
ensembles = randn(N*N, num_ensembles); % Initial ensembles
% Ensemble mean
mean_ensemble = mean(ensembles, 2);
m(:,1) = mean_ensemble;

timer = tic;

for i = 1:n_t
    % Forecast step: propagate each ensemble member
    for k = 1:num_ensembles
        ensemble_k = ensembles(:,k);
        % Propagate using explicit Euler
        ensemble_k = ensemble_k + dt*(-Gamma*ensemble_k + Sigma*randn(N*N,1));
        ensembles(:,k) = ensemble_k;
    end
    
    % Analysis step: update with observations
    if mod(i, nsavef) == 0
        % Observation matrix
        H = zeros(2*Ltracers, N*N);
        for n = 1:Ltracers
            Xh = exp(1i*(kxx*X((n-1)*2 + 1, i) + kyy*X((n-1)*2 + 2, i)));
            x1 = Xh.*projx;
            x2 = Xh.*projy;
            H((n-1)*2 + 1,:) = x1(:).';
            H((n-1)*2 + 2,:) = x2(:).';
        end
        d = dXdt(:,i);
        y = H * ensembles; % Observations of ensembles
        
        % Calculate Kalman gain
        P_yy = (y * y') / (num_ensembles - 1) + sigmax^2 * eye(2*Ltracers);
        P_xy = (ensembles * y') / (num_ensembles - 1);
        K = P_xy / P_yy;
        
        % Update ensembles
        for k = 1:num_ensembles
            ensembles(:,k) = ensembles(:,k) + K * (d - y(:,k));
        end
        
        % Save ensemble mean
        mean_ensemble = mean(ensembles, 2);
        m(:,i/nsavef+1) = mean_ensemble;
    end
end

fprintf('filter time %s\n', time2str(toc(timer)));

%-------------------------------------------------------------------------%
% save only dt_save time steps

Xt = X(:,1:nsavef:length(t)); 
pin.computed = pout; % save computed parameters

casestr = casename(pin.N, pin.Ltracers, sigmax, pin.Tend, pin.spectrum_type, pin.damping_type, pin.dtsave);
datafolder = fullfile(datafolder,casestr);
if fsavedata && ~isfolder(datafolder); mkdir(datafolder); end

if fsavedata
    save(fullfile(datafolder, 'data.mat'), 'm', 'Xt', 'pin', '-v7.3')
end

end
