% Clear workspace
clear;

% Load cases and select one
cases = loadcases();
casename = cases(3);
datafolder = fullfile('data', casename);

% Load external data
dataext = load(fullfile(datafolder, 'diag11.mat'));

% Load time data
fnamet = fullfile(datafolder, 't11.bin'); 
fidt = fopen(fnamet, 'r');
t = fread(fidt, 'double'); 
fclose(fidt);
dt = t(2) - t(1);
nt = length(t);
Ntruth = dataext.params.N;

% Load observations
Ltracers = 10;
sigmax = 0.1;
obsdatafile = dir(fullfile(sprintf('obs_Ltracers%i_sigmax%.2g.mat', Ltracers, sigmax)));
if isempty(obsdatafile)
    error('obs data not found');
end
obsdata = load(fullfile(obsdatafile.folder, obsdatafile.name));

% Set up the pin structure with necessary parameters
pin.Tend = t(end);
pin.Ltracers = Ltracers;
pin.N = 8; % Grid points in each dimension (must be even)
pin.X = obsdata.X;
pin.E0 = 0.8;
pin.damping_type = 'unif';
pin.d = 0.1;
pin.mu = 0.05;
pin.spectrum_type = 'kolmo';
pin.num_modes_flat = 2;
pin.sigmax = obsdata.sigmax;
pin.dt = dt;
pin.dtsave = dt * 10;
pin.fsavedata = true;
pin.datafolder = 'Ahmad_QG';

pin.beta = 0.05;
pin.dd = 5E-3;
pin.epsf = 0;
pin.kf = 50;  % No effect since epsf = 0
pin.sigf = 1; % No effect since epsf = 0
pin.alpf = 0; % No effect since epsf = 0
pin.tau0 = 0.04;
pin.tdel = 0; % No effect for sin forcing
pin.forcingtype = 'sin';

Nt = 2e5;

% Run the Particle Filter with QG Dynamics
[pin, m, Xt] = LDA_PF_QG_withDynamics(pin);

% Coarsen qh data
timer = tic;
qhc = coarsenqh(datafolder, pin.dtsave);
fprintf('coarsenqh time %s\n', time2str(toc(timer)));

% Set plot options
opts.renderer = 'zbuffer';
opts.plot_times = [100 200 300];
opts.fsavefigs = true;

% Perform analysis
analysis(pin, opts, m, qhc, Xt);

% Particle Filter with QG Dynamics function
function [pin, m, Xt] = LDA_PF_QG_withDynamics(pin)
    % Extract input parameters
    Tend = pin.Tend;
    damping_type = pin.damping_type;  
    d = pin.d;
    mu = pin.mu;
    spectrum_type = pin.spectrum_type;
    E0 = pin.E0;
    num_modes_flat = 2;
    sigmax = pin.sigmax;
    Ltracers = pin.Ltracers;
    N = pin.N; 
    dt = pin.dt;
    dtsave = pin.dtsave;
    fsavedata = pin.fsavedata; 
    datafolder = pin.datafolder;
    nu = pin.nu; % Coefficient of biharmonic vorticity diffusion
    sh = pin.sh; % Hyperdiffusion exponent
    kf = pin.kf; % Forcing wavenumber
    epsf = pin.epsf; % Forcing amplitude
    sigf = pin.sigf; % Forcing width
    alpf = pin.alpf; % Forcing anisotropy parameter

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
    weights = ones(ensemble_size, 1) / ensemble_size; % Initial weights
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

        % QG Dynamics Forecast for the ensemble
        for j = 1:ensemble_size
            q = reshape(X_ens(:, j), [N, N]);
            
            % Update q using QG dynamics
            q = QG_Dynamics(q, N, dt, nu, sh, kf, epsf, sigf, alpf); % Replace with actual QG step function
            
            % Reshape q back to vector form
            X_ens(:, j) = q(:);
        end

        % Particle Filter update step
        obs = dXdt(:, i); % Observations at the current time step
        for j = 1:ensemble_size
            pred_obs = Px * X_ens(:, j); % Predicted observation
            weights(j) = weights(j) * exp(-0.5 * (obs - pred_obs)' * (obs - pred_obs) / sigmax^2);
        end
        weights = weights / sum(weights); % Normalize weights

        % Resampling step
        resample_indices = randsample(ensemble_size, ensemble_size, true, weights);
        X_ens = X_ens(:, resample_indices);
        weights = ones(ensemble_size, 1) / ensemble_size; % Reset weights after resampling

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
        save(fullfile(datafolder, 'qg_data.mat'), 'm', 'Xt', 'pin', '-v7.3')
    end
end

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
    end
end

% Function for QG Dynamics
function q = QG_Dynamics(q, N, dt, nu, sh, kf, epsf, sigf, alpf)
    % Setup hyperviscous PV dissipation and linear damping
    k = [0:N/2 -N/2+1:-1]'; % wavenumbers
    L = zeros(N, N);
    for jj = 1:N
        for ii = 1:N
            kr = sqrt(k(ii)^2 + k(jj)^2);
            L(ii,jj) = -nu * kr^(2 * sh);
        end
    end
    
    % Placeholder for RHS_Spectral function
    k0 = RHS_Spectral(q, N, nu, sh, kf, epsf, sigf, alpf);
    l0 = L .* q;
    
    % ARK4 Stages
    M = 1 ./ (1 - 0.25 * dt * L);
    
    q1 = M .* (q + 0.5 * dt * k0 + 0.25 * dt * l0);
    k1 = RHS_Spectral(q1, N, nu, sh, kf, epsf, sigf, alpf);
    l1 = L .* q1;
    
    q2 = M .* (q + dt * (13861 * k0 / 62500 + 6889 * k1 / 62500 + 8611 * l0 / 62500 - 1743 * l1 / 31250));
    k2 = RHS_Spectral(q2, N, nu, sh, kf, epsf, sigf, alpf);
    l2 = L .* q2;
    
    q3 = M .* (q + dt * (-0.04884659515311858 * k0 - 0.1777206523264010 * k1 + 0.8465672474795196 * k2 + 0.1446368660269822 * l0 - 0.2239319076133447 * l1 + 0.4492950415863626 * l2));
    k3 = RHS_Spectral(q3, N, nu, sh, kf, epsf, sigf, alpf);
    l3 = L .* q3;
    
    q4 = M .* (q + dt * (-0.1554168584249155 * k0 - 0.3567050098221991 * k1 + 1.058725879868443 * k2 + 0.3033959883786719 * k3 + 0.09825878328356477 * l0 - 0.5915442428196704 * l1 + 0.8101210538282996 * l2 + 0.2831644057078060 * l3));
    k4 = RHS_Spectral(q4, N, nu, sh, kf, epsf, sigf, alpf);
    l4 = L .* q4;
    
    q5 = M .* (q + dt * (0.2014243506726763 * k0 + 0.008742057842904184 * k1 + 0.1599399570716811 * k2 + 0.4038290605220775 * k3 + 0.2260645738906608 * k4 + 0.1579162951616714 * l0 + 0.1867589405240008 * l2 + 0.6805652953093346 * l3 - 0.2752405309950067 * l4));
    k5 = RHS_Spectral(q5, N, nu, sh, kf, epsf, sigf, alpf);
    l5 = L .* q5;

    % Stochastic forcing
    dW1 = randn(N, N/2+1);
    dW2 = randn(N, N/2+1);
    dW = (dW1 + 1i * dW2) / sqrt(2);
    dW(1,1) = 0;
    dW(1, N/2+1) = dW1(1, N/2+1);
    dW(N/2+1, 1) = dW1(N/2+1, 1);
    dW(N/2+1, N/2+1) = dW1(N/2+1, N/2+1);
    dW(N/2+2:end, 1) = conj(dW(N/2:-1:2, 1));
    dW(N/2+2:end, N/2+1) = conj(dW(N/2:-1:2, N/2+1));
    
    Fk = compute_Fk(N, epsf, sigf, alpf, kf); % Compute Fk based on parameters

    xik = N^2 * [dW, conj(dW([1, N:-1:2], N/2:-1:2))] .* sqrt(Fk);

    % Successful step, proceed to evaluation
    qp = real(ifft2(q + dt * (0.1579162951616714 * (k0 + l0) + 0.1867589405240008 * (k2 + l2) + 0.6805652953093346 * (k3 + l3) - 0.2752405309950067 * (k4 + l4) + (k5 + l5) / 4) + sqrt(dt) * xik));
    q = fft2(qp);
end


function Fk = compute_Fk(N, epsf, sigf, alpf, kf)
    % Compute spectral forcing Fk
    Fk = zeros(N, N);
    k = [0:N/2 -N/2+1:-1]'; % wavenumbers
    for ii = 1:N
        for jj = 1:N
            kr = sqrt(k(ii)^2 + k(jj)^2);
            phi = atan(k(jj) / k(ii));
            if isnan(phi)
                phi = 0;
            end
            ampl = epsf * (1 + alpf * cos(2 * phi));
            Fk(ii, jj) = (ampl * kr) * exp(-0.5 * (kr - kf)^2 / sigf^2) / (sqrt(2 * pi) * sigf);
        end
    end
    Fk(1, 1) = 0;
end

% Helper function for initial name generation
function casestr = casename(N, Ltracers, sigmax, Tend, spectrum_type, damping_type, dtsave)
    casestr = sprintf('N%i_L%i_sigmax%.2f_T%.1f_spectrum%s_damping%s_dtsave%.2f', ...
                      N, Ltracers, sigmax, Tend, spectrum_type, damping_type, dtsave);
end

