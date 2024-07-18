function [pin, m, Xt] = Algorithm_design(pin)
% Ensuring all necessary fields are initialized
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fields = {'Tend', 'damping_type', 'd', 'mu', 'spectrum_type', 'E0', ...
          'num_modes_flat', 'Ltracers', 'N', 'dt', 'dtsave', 'fsavedata', 'datafolder', 'X'};
defaults = {10, 'unif', 0.1, 0, 'kolmo', 1, 2, 10, 32, 0.01, 0.1, false, '', randn(32*32, 1)};

for i = 1:length(fields)
    if ~isfield(pin, fields{i})
        pin.(fields{i}) = defaults{i};
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define parameters from pin for clarity
[Tend, damping_type, d, mu, spectrum_type, E0, num_modes_flat, Ltracers, N, ...
 dt, dtsave, fsavedata, datafolder, X] = deal(pin.Tend, pin.damping_type, pin.d, pin.mu, ...
 pin.spectrum_type, pin.E0, pin.num_modes_flat, pin.Ltracers, pin.N, pin.dt, ...
 pin.dtsave, pin.fsavedata, pin.datafolder, pin.X);

% Time and grid setup
t = 0:dt:Tend;
n_t = length(t);
dtsave_indices = 1:round(dtsave/dt):n_t;
tsave = t(dtsave_indices);

% Grid and spectral setup
L = 2*pi;
x = linspace(0, L, N+1);
x(end) = []; % periodic boundary condition
y = x;
[xx, yy] = meshgrid(x, y);
kx = fftshift(-N/2:N/2-1);
ky = kx;
[kxx, kyy] = meshgrid(kx, ky);
knorm = hypot(kxx, kyy);
knorm(knorm == 0) = 1; % Avoid division by zero for zero mode

% Damping and spectral setup
gamma = d + mu*knorm.^2; % Example damping


N = pin.N;
n_ens = 3;  % Number of ensemble members
sigma_noise = 0.01;  % Noise level for the forecast model

% Initialize ensemble members by perturbing pin.X
X_ens = cell2mat(pin.X) ;

% Define Sigma and Gamma for the Kalman filter update step
Sigma = 0.1 * eye(size(pin.X, 1));  % Measurement noise covariance
Gamma = 0.01 * eye(size(pin.X, 1));  % Process noise covariance

for k = 1:n_t

    display(k);

    for j = 1:n_ens
        X_ens(:, j) = forecast_model(X_ens(:, j), dt, N, L);
    end

    observation = pin.X ;
    
    X_ens = update_step(X_ens, observation, Sigma, Gamma, N);

    m(:, k) = mean(X_ens, 2);

    if any(k == dtsave_indices)
        Xt(:, k == dtsave_indices) = m(:, k);

    end
end

if fsavedata
    casestr = casename(N, Ltracers, sigma_noise, Tend, spectrum_type, damping_type, dtsave);
    datafolder_full = fullfile(datafolder, casestr);
    if ~exist(datafolder_full, 'dir')
        mkdir(datafolder_full);
    end
    save(fullfile(datafolder_full, 'hi_data.mat'), 'm', 'Xt', 'pin');
end

function X_next = forecast_model(X_current, dt, N, L)

    A = diag(0.95 * ones(N*N, 1));  % Linear damping matrix
    noise = sigma_noise * randn(N*N, 1);  % Additive Gaussian noise
    X_next = X_current;   % Simple Euler integration
end

function X_ens = update_step(X_ens, observation, Sigma, Gamma, N)

    P_forecast = cov(X_ens');  % Forecast covariance
    K = P_forecast * (P_forecast )^(-1);  % Kalman gain
    observation_expanded = repmat(observation, 1, size(X_ens, 2));
    innovations = observation_expanded;
    X_ens = X_ens;
end

end