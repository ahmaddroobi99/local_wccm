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

% QG parameters
beta = pin.beta;
kf = pin.kf;
epsf = pin.epsf;
sigf = pin.sigf;
alpf = pin.alpf;
tau0 = pin.tau0;
tdel = pin.tdel;
forcingtype = pin.forcingtype;
dd = pin.dd;

% Set physical parameters
sh = 8; % hyperdiffusion exponent (needs to be even)
if N == 128
    if epsf == 0
        nu = 2.5e-28 *(50/50)^(-2*sh+2/3); % Coefficient of biharmonic vorticity diffusion
    else
        nu = 2.5e-28 *(kf/50)^(-2*sh+2/3); % Coefficient of biharmonic vorticity diffusion
    end
elseif N == 256
    nu = 5E-33 *(kf/50)^(-2*sh+2/3);
    if tau0 == 0
        nu = nu*100;
    end
elseif N == 512
    nu = 1E-34 *(kf/50)^(-2*sh+2/3);
else
    % Add a default value for nu if N is not one of the predefined values
    nu = 1E-34 *(kf/50)^(-2*sh+2/3);
end

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
% damping_type = ' ' ;
% spectrum_type = ' ' ;  


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

% Set up hyperviscous PV dissipation and linear damping
k = [0:N/2 -N/2+1:-1]'; % wavenumbers
L = zeros(N,N);
for jj=1:N
    for ii=1:N
        kr = sqrt(k(ii)^2+k(jj)^2);
        L(ii,jj) = -nu*kr^(2*sh);
    end
end

% Set up stochastic forcing with spectral
Fk = zeros(N,N);
for ii=1:N
    for jj=1:N
        kr = sqrt(k(ii)^2+k(jj)^2);
        phi = atan(k(jj)/k(ii));
        if isnan(phi)
            phi = 0;
        end
        ampl = epsf*(1+alpf*cos(2*phi));
        Fk(ii,jj) = (ampl*kr)*exp(-.5*(kr-kf)^2/sigf^2)/(sqrt(2*pi)*sigf);
    end
end
Fk(1,1) = 0;

% Surface wind stress
dx = 2*pi/N;
[KX, KY] = meshgrid(k,k);
[X_mesh, Y_mesh] = meshgrid(-pi:dx:pi-dx,-pi:dx:pi-dx);

if strcmp(forcingtype,'sin')
    tau = tau0 *sin(Y_mesh).^2;
elseif strcmp(forcingtype,'sinh')
    tau = tau0 *sech(Y_mesh/tdel).^2;
elseif strcmp(forcingtype,'constant')
    tau = tau0;
end
tauk = -1i*KY.*fft2(tau);
clear k KX KY

% Initialize diagnostics
countDiag = 100; % Compute diagnostics every countDiag steps
num_diag_points = floor(n_t / countDiag); % Ensure this is an integer
T = zeros(1, num_diag_points);
vb = zeros(1, num_diag_points);
utz = zeros(N, num_diag_points);
energy = zeros(N/2+1, num_diag_points);
enstrophy = zeros(N/2+1, num_diag_points);

qk1 = zeros(N,N); % mean
qk2 = zeros(N,N); % var
count = 0;

timer = tic;             % Start timer
Px = zeros(2 * Ltracers, N * N); % Matrix for projections
Sigma2 = Sigma * Sigma'; % Covariance of the model error

% Main Loop
for i = 1:n_t
    % Compute projection matrix Px
    for n = 1:Ltracers
        if (n-1)*2 + 2 <= size(X, 1) && i <= size(X, 2)
            fprintf('Computing projection for tracer %d at time step %d\n', n, i); % Debug statement
            Xh = exp(1i * (kxx * X((n-1)*2 + 1, i) + kyy * X((n-1)*2 + 2, i)));
            Px((n-1)*2 + 1, :) = reshape(Xh .* projx, 1, []);
            Px((n-1)*2 + 2, :) = reshape(Xh .* projy, 1, []);
        else
            fprintf('Index exceeds matrix dimensions at tracer %d and time step %d\n', n, i); % Debug statement
            continue; % Skip the iteration if indices are out of bounds
        end
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

        % Add QG dynamics update step here
        q = QG_Step(q, nu, kx, ky, dt, tauk, L, beta); % Ensure 'tauk' is correctly computed

        X_ens(:, j) = q(:);
    end

    % Observation update for each ensemble member
    for j = 1:ensemble_size
        X_ens(:, j) = X_ens(:, j) + K * (dXdt(:, i) - Px * X_ens(:, j));
    end

    % Save results at specified intervals
    if mod(i, nsavef) == 0
        m(:, i/nsavef + 1) = X_mean; % Store the mean of the ensemble
    end
end

fprintf('Time for simulation: %f seconds\n', toc(timer));
pin.m = m;
Xt = X; % Return observed data as Xt (for compatibility)
end

% Supporting Functions

function gamma = get_gamma(damping_type, d, mu, knorm)
    if strcmp(damping_type, 'constant')
        gamma = d * ones(size(knorm));
    elseif strcmp(damping_type, 'unif')
        gamma = d * (knorm / max(knorm));
    elseif strcmp(damping_type, 'viscous')
        gamma = d * (knorm.^mu / max(knorm));
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
    end
end

function idx = zeromodes2d(N)
    idx = false(N, N);
    idx(1, 1) = true;
end

function idx = junkmodes2d(N)
    idx = false(N, N);
    for ii = 1:N
        for jj = 1:N
            if hypot(ii - N/2 - 1, jj - N/2 - 1) > N/2
                idx(ii, jj) = true;
            end
        end
    end
end

function [energy, enstrophy] = compute_energy_enstrophy(qk, kxx, kyy)
    k2 = kxx.^2 + kyy.^2;
    energy = sum(abs(qk).^2 .* k2, 'all') / 2;
    enstrophy = sum(abs(qk).^2 .* k2.^2, 'all') / 2;
end

% Placeholder for the QG dynamics update function
function q = QG_Step(q, nu, kx, ky, dt, tauk, L, beta)
    % Example QG dynamics update step
    % Here, implement the specific QG dynamics update for the model state q

    % Fourier transform of q
    q_hat = fft2(q);

    % Update q_hat using the QG model dynamics
    % Placeholder for the QG model dynamics update
    % Implement the actual QG dynamics here
    q_hat_new = q_hat; % Modify this line with the actual QG update

    % Inverse Fourier transform to obtain the updated q
    q = real(ifft2(q_hat_new));
end
