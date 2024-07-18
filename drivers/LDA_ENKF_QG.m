% Key Points:

% Initialization: Initialize the ensemble of model states and set up the grid and wave vectors for the QG model.

% Observation Model: The observations of the Lagrangian tracers are used to update the model state.

% Ensemble Prediction: Each ensemble member is evolved using the QG dynamics. Here, you need to implement the specific update steps for the QG equations.

% EnKF Update: The ensemble is updated by assimilating the observed tracer positions using the Kalman gain.



% Assumptions and Modifications:

% The QG_Step function or equivalent update mechanism for the QG model dynamics needs to be defined based on  specific model setup (under Research).
% Adjust the number of tracers, ensemble size, and other parameters according to specific problem requirements.
%  proper initialization and observation handling for specific case.

% ||||||||||||||||||||||||||||||||||||||||||||||||||||
% Implementation





function [pin, m, Xt] = LDA_ENKF_QG(pin)
% LDA_ENKF_QG: Implementation of the Ensemble Kalman Filter (EnKF) for Lagrangian Data Assimilation in the Quasi-Geostrophic (QG) model.

% Extract input parameters
Tend = pin.Tend;
damping_type = pin.damping_type;  
d = pin.d;
mu = pin.mu;
spectrum_type = pin.spectrum_type;
E0 = pin.E0;
num_modes_flat = 2;
sigmax = 0.1; % Assumed constant value, can be parameterized as needed
Ltracers = pin.Ltracers;
N = pin.N; 
dt = pin.dt;
dtsave = pin.dtsave;
fsavedata = pin.fsavedata; 
datafolder = pin.datafolder;

% Setup simulation parameters
t = 0:dt:Tend+dt;
n_t = length(t);
dtsave = round(dtsave/dt)*dt; % Round to a multiple of dt
nsavef = dtsave/dt;
tsave = t(1:nsavef:length(t));
n_tsave = length(tsave);

fprintf('Total modes (2*N+1) = %i^2, max wavenumber resolved %d\n', N, N/2-1);
Ntotal = N*N;

% Initialize grids and wave vectors
L = 2*pi;
x = L/N*(0:N-1)';
y = L/N*(0:N-1)';
[xx, yy] = meshgrid(x, y);
kx = [0:N/2-1 0 -N/2+1:-1];
ky = [0:N/2-1 0 -N/2+1:-1];
[kxx, kyy] = meshgrid(kx, ky);
knorm = hypot(kxx, kyy);

projx = (1i * -kyy ./ knorm);
projy = (1i * kxx ./ knorm);
projx(zeromodes2d(N)) = 0;
projy(zeromodes2d(N)) = 0;

% Damping and spectrum types
gamma = get_gamma(damping_type, d, mu, knorm);
sigma = get_sigma(spectrum_type, E0, knorm, gamma, num_modes_flat, N, sigmax);

sigma(junkmodes2d(N)) = 0;
Sigma = sigma(:) .* speye(Ntotal);
gamma(junkmodes2d(N)) = 0;
Gamma = gamma(:).' .* speye(Ntotal);

% Save computed parameters
pout = struct('tsave', tsave, 'dtsave', dtsave, 'N', N, 'Ntotal', Ntotal, ...
    'sigma', sigma, 'Sigma', Sigma, 'gamma', gamma, 'Gamma', Gamma, ...
    'x', x, 'y', y, 'kx', kx, 'ky', ky, 'kxx', kxx, 'kyy', kyy);
pin.computed = pout;

% Initial ensemble setup
ensemble_size = 250;
mprev = zeros(N*N, 1);
X_ens = mvnrnd(mprev', eye(N*N), ensemble_size)';
m = zeros(N*N, n_tsave); % Array to store the results
X = cell2mat(pin.X);     % Observation matrix
dXdt = [(X(:,2)-X(:,1)) (X(:,2:end) - X(:,1:end-1))] / dt; % Time derivative of observations

timer = tic;             % Start timer
Px = zeros(2 * Ltracers, N * N); % Matrix for projections
Sigma2 = Sigma * Sigma'; % Covariance of the model error

% Main Loop
for i = 1:n_t
    % Compute projection matrix Px
    for n = 1:Ltracers
        Xh = exp(1i * (kxx * X((n-1)*2 + 1, i) + kyy * X((n-1)*2 + 2, i)));
        Px((n-1)*2 + 1, :) = reshape(Xh .* projx, 1, []);
        Px((n-1)*2 + 2, :) = reshape(Xh .* projy, 1, []);
    end

    % Ensemble Kalman Filter (EnKF) update
    X_mean = mean(X_ens, 2);  % Mean of the ensemble
    X_prime = X_ens - X_mean; % Deviations from the mean
    C = X_prime * X_prime' / (ensemble_size - 1); % Ensemble covariance

    % Compute the Kalman gain
    Cov = Px * C * Px' + sigmax^-2 * eye(size(Px, 1)); 
    K = C * Px' / Cov;

    % QG Dynamics Forecast for the ensemble
    for j = 1:ensemble_size
        q = reshape(X_ens(:, j), [N, N]);
        
        % Update q using QG dynamics
        % Add your QG dynamics update step here
        % Example: q = QG_Step(q, params); 
        % Replace QG_Step with the actual function or code snippet that updates q
         
        % Here, `QG_Step` is a placeholder for your actual QG dynamics update code
        % Example: q = q + dt * (nonlinear_term + diffusion_term + external_forcing);
        
        % Reshape q back to vector form
        X_ens(:, j) = q(:);
    end

    % Apply the EnKF update
    X_ens_diff = dXdt(:, i) - Px * X_ens;
    X_ens_update = X_ens + K * X_ens_diff;

    % Update the EnKF ensemble
    X_ens = X_ens_update;
    X_ens_reshaped = mean(X_ens, 2);

    if mod(i, nsavef) == 0
        m(:, floor(i/nsavef) + 1) = X_ens_reshaped;
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
% % % % % % % % % % % % % % % % % % % % % ----------------------------------------------
% Helper function to get gamma based on damping type
function gamma = get_gamma(damping_type, d, mu, knorm)
    if strcmp(damping_type, 'scaledep')
        gamma = d + mu * knorm.^2;
    elseif strcmp(damping_type, 'unif')
        gamma = d * ones(size(knorm));
    else
        error('Unknown damping type');
    end
end

function sigma = get_sigma(spectrum_type, E0, knorm, gamma, num_modes_flat, N, sigmax)
    if strcmp(spectrum_type, 'kolmo')
        sigma = sqrt(2 * gamma .* (E0 * knorm.^(-5/3)) ./ knorm);
    elseif strcmp(spectrum_type, 'com_kolmo')
        sigma = sqrt(2 * gamma .* (E0 * (knorm / num_modes_flat).^(-5/3)) ./ knorm);
    elseif strcmp(spectrum_type, 'unif')
        sigma = sigmax * ones(size(knorm));
    elseif strcmp(spectrum_type, 'none')
        sigma = zeros(size(knorm));
    else
        error('Unknown spectrum type');
% _________________________



    end
end
