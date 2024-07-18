function [pin, m, Xt] = runQGCGEnkf(pin)
% inputs
Tend = pin.Tend;

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







% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Solving the filter for velocity recovery (EXPLICIT EULER)
% Embedded with Ensemble Kalman Filter (EnKF)

% Initialization
is = 1;
ensemble_size = 300;
Rprev = zeros(N*N, N*N); % Initial covariance matrix
mprev = zeros(N*N, 1);    % Initial state vector
X_ens = mvnrnd(mprev', Rprev, ensemble_size)'; % Initial ensemble

m = zeros(N*N, n_tsave); % Array to store the results
X = cell2mat(pin.X);     % Observation matrix

dXdt = [(X(:,2)-X(:,1)) (X(:,2:end) - X(:,1:end-1))] / dt; % Time derivative of observations

timer = tic;             % Start timer
Px = zeros(2*Ltracers, N*N); % Matrix for projections
Sigma2 = Sigma * Sigma'; % Covariance of the model error

% Main Loop
for i = 1:n_t
    % Compute projection matrix Px
    for n = 1:Ltracers
        Xh = exp(1i * (kxx * X((n-1)*2 + 1, i) + kyy * X((n-1)*2 + 2, i)));
        x1 = Xh .* projx;
        x2 = Xh .* projy;
        Px((n-1)*2 + 1, :) = x1(:);
        Px((n-1)*2 + 2, :) = x2(:);
    end

    % Gaussian nonlinear filter equations
    RpPxt = Rprev * Px';
    GammaRprev = Gamma * Rprev;
    Rtemp = Rprev + dt * (-GammaRprev - GammaRprev' + Sigma2 - sigmax^-2 * (RpPxt * RpPxt'));
    mtemp = mprev + dt * ((-Gamma * mprev + F) + sigmax^-2 * RpPxt * (dXdt(:, i) - Px * mprev));

    % Ensemble Kalman Filter (EnKF) update
    X_mean = mean(X_ens, 2);  % Mean of the ensemble
    X_prime = X_ens - X_mean; % Deviations from the mean
    C = X_prime * X_prime' / (ensemble_size - 1); % Ensemble covariance

    % Compute the Kalman gain
    Cov = Px * C * Px' + sigmax^-2 * eye(size(Px, 1)); 
    K = RpPxt / Cov;

    % Forecast the ensemble
    X_ens_forecast = X_ens + dt * (-Gamma * X_ens + F) + sqrt(dt) * randn(size(X_ens));

    % Apply the EnKF update
    X_ens_diff = dXdt(:, i) - Px * X_ens_forecast;
    X_ens_update = X_ens_forecast + K * X_ens_diff;

    % Update state estimates
    mprev = mtemp;
    Rprev = Rtemp;

    % Update the EnKF ensemble
    X_ens = X_ens_update;
    X_ens_reshaped = mean(X_ens, 2);

    if mod(i, nsavef) == 0
        is = is + 1;
        m(:, is) = X_ens_reshaped;
    end
end

fprintf('Filter time %s\n', time2str(toc(timer)))

% Save only dt_save time steps
Xt = X(:, 1:nsavef:length(t)); 

pin.computed = pout; % Save computed parameters

casestr = casename(pin.N, pin.Ltracers, sigmax, pin.Tend, pin.spectrum_type, pin.damping_type, pin.dtsave);
datafolder = fullfile(datafolder, casestr);

if fsavedata && ~isfolder(datafolder)
    mkdir(datafolder);
end

if fsavedata
    save(fullfile(datafolder, 'hello_data.mat'), 'm', 'Xt', 'pin', '-v7.3')
end


end



