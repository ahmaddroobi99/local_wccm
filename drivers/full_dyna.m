function [pin, m, Xt] = LDA_ENKF_QG(pin)
    % LDA_ENKF_QG: Implementation of the Ensemble Kalman Filter (EnKF)
    % for Lagrangian Data Assimilation in the Quasi-Geostrophic (QG) model.

    % Extract input parameters
    Tend = pin.Tend;
    damping_type = pin.damping_type;  
    d = pin.d;
    mu = pin.mu;
    spectrum_type = pin.spectrum_type;
    E0 = pin.E0;   
    num_modes_flat = 2;
    sigmax = 0.1; % Standard deviation of observation noise
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

    % QG parameters
    qlim = 1.5E5; % if any q > qlim, simulation stops
    beta = 0.05;
    dd = 5E-3;

    epsf = 0;
    kf = 50;  % no effect since epsf = 0
    sigf = 1; % no effect since epsf = 0
    alpf = 0; % no effect since epsf = 0
    tau0 = 0.04;
    tdel = 0; % no effect for sin forcing
    forcing_type = 'sin';

    Nt = 2e5;
    % Set physical parameters
    sh = 8; % hyperdiffusion exponent (needs to be even)
    if N == 8
        if epsf == 0
            nu = 2.5e-28 *(50/50)^(-2*sh+2/3); %5E-33; % Coefficient of biharmonic vorticity diffusion
        else
            nu = 2.5e-28 *(kf/50)^(-2*sh+2/3); %5E-33; % Coefficient of biharmonic vorticity diffusion
        end
    elseif N == 256
        nu = 5E-33 *(kf/50)^(-2*sh+2/3); %5E-33; % Coefficient of biharmonic vorticity diffusion
        if tau0 == 0
            nu = nu*100;
        end
    elseif N == 512
        nu = 1E-34 *(kf/50)^(-2*sh+2/3);
    end

    % Put useful stuff into a struct    
    params = struct('beta', beta, 'sh', sh, 'dd', dd,'N', N, 'nu', nu, ...
                    'kf', kf, 'epsf', epsf, 'sigf', sigf, 'alpf', alpf, ...
                    'tau0', tau0, 'tdel', tdel, 'forcingtype', forcing_type);

    F = 0; % zero forcing or periodic forcing only (not tested)

    % Setup grid and wave vectors
    L = 2*pi;
    x = L/N*(0:N-1)'; % x coordinate
    y = L/N*(0:N-1)'; % y coordinate
    [xx, yy] = meshgrid(x, y); % xy grid

    kx = [0:N/2-1 0 -N/2+1:-1]; % x wave vector
    ky = [0:N/2-1 0 -N/2+1:-1]; % y wave vector
    [kxx, kyy] = meshgrid(kx, ky);
    knorm = hypot(kxx, kyy);

    projx = (1i * -kyy ./ knorm);
    projy = (1i * kxx ./ knorm);
    projx(zeromodes2d(N)) = 0;
    projy(zeromodes2d(N)) = 0;

    % Determine gamma based on damping type
    if strcmp(damping_type, 'scaledep')
        gamma = (d + mu * knorm.^2);
    elseif strcmp(damping_type, 'unif')
        gamma = d * ones(N, N);
    else
        error('Unknown damping type')
    end

    % Determine sigma based on spectrum type
    if strcmp(spectrum_type, 'kolmo')
        sigma = sqrt(2 * gamma .* (E0 * knorm.^(-5/3)) ./ knorm);
    elseif strcmp(spectrum_type, 'com_kolmo')
        sigma = sqrt(2 * gamma .* (E0 * (knorm / num_modes_flat).^(-5/3)));
        sigma(knorm < num_modes_flat) = sqrt(2 * gamma(knorm < num_modes_flat) * E0 ./ knorm(knorm < num_modes_flat));
    elseif strcmp(spectrum_type, 'third')
        sigma = sqrt(2 * gamma .* (E0 * knorm.^(-3)) ./ knorm);
    elseif strcmp(spectrum_type, 'com_third')
        sigma = sqrt(2 * gamma .* (E0 * (knorm / num_modes_flat).^(-3)));
        sigma(knorm < num_modes_flat) = sqrt(2 * gamma(knorm < num_modes_flat) * E0 ./ knorm(knorm < num_modes_flat));
    elseif strcmp(spectrum_type, 'fourth')
        sigma = sqrt(2 * gamma .* (E0 * knorm.^(-4)) ./ knorm);
    elseif strcmp(spectrum_type, 'flat')
        sigma = sqrt(2 * gamma .* (E0 * knorm.^(0)) ./ knorm);
    elseif strcmp(spectrum_type, 'exp')
        if isscalar(pin.sigma)
            sigma = pin.sigma * ones(N, N);
        else
            sigma = pin.sigma;
        end
    else
        error('Unknown spectrum type')
    end

    sigma(junkmodes2d(N)) = 0;
    Sigma = sigma(:) .* speye(Ntotal);
    gamma(junkmodes2d(N)) = 0;
    Gamma = gamma(:) .* speye(Ntotal);

    % Save computed parameters
    pout = struct('tsave', tsave, 'dtsave', dtsave, 'N', N, 'Ntotal', Ntotal, ...
                  'sigma', sigma, 'Sigma', Sigma, 'gamma', gamma, 'Gamma', Gamma, ...
                  'x', x, 'y', y, 'kx', kx, 'ky', ky, 'kxx', kxx, 'kyy', kyy);
    pin.computed = pout;
%-------------------------------------------------------------------------%
% Solve filter for velocity recovery using QG dynamics (EXPLICIT EULER)
%-------------------------------------------------------------------------%

% Initial setup
m = zeros(N*N, n_tsave);
X = cell2mat(pin.X);
dXdt = [(X(:,2)-X(:,1)) (X(:,2:end) - X(:,1:end-1))]./dt;

% Initial conditions
is = 1;
Rprev = zeros(N*N, N*N); %1/2*(Sigma*Sigma')./Gamma;
mprev = zeros(N*N,1);
m(:,is) = mprev;

timer = tic;
Px = zeros(2*Ltracers, N*N);

Sigma2 = Sigma*Sigma';

% Main loop
for i = 1:n_t
    % Compute projection matrix Px
    for n = 1:Ltracers
        Xh = exp(1i*(kxx*X((n-1)*2 + 1, i) + kyy*X((n-1)*2 + 2, i)));
        x1 = Xh.*projx;
        x2 = Xh.*projy;
        Px((n-1)*2 + 1,:) = x1(:).';
        Px((n-1)*2 + 2,:) = x2(:).';
    end     
    RpPxt = Rprev*Px';
    GammaRprev = Gamma*Rprev;

    % Update mean with QG dynamics
    mtemp = mprev + dt*((-Gamma*mprev + F) + sigmax^-2*RpPxt*(dXdt(:,i) - Px*mprev));
    mtemp_reshaped = reshape(mtemp, N, N);
    mtemp_qg = QG_dynamics(mtemp_reshaped, dt, params);
    mtemp = reshape(mtemp_qg, [], 1);

    % Update mprev
    mprev = mtemp;

    if mod(i, nsavef) == 0
        is = is + 1;
        m(:,is) = mprev;
    end
end

fprintf('Filter time %s\n', time2str(toc(timer)));



% Function to compute QG dynamics
function q_next = QG_dynamics(q, dt, params)
    N = params.N;
    tauk = compute_tauk(params, size(q));

    % Set up hyperviscous PV dissipation and linear damping
    k = [0:N/2 -N/2+1:-1]'; % wavenumbers
    LL = zeros(N, N);
    for jj = 1:N
        for ii = 1:N
            kr = sqrt(k(ii)^2 + k(jj)^2);
            LL(ii, jj) = -params.nu * kr^(2 * params.sh);
        end
    end

    % Initialize
    M = 1 ./ (1 - .25 * dt * LL);

    % First stage ARK4
    k0 = RHS_Spectral(q, params) + tauk;
    l0 = LL .* q;

    % Second stage
    q1 = M .* (q + .5 * dt * k0 + .25 * dt * l0);
    k1 = RHS_Spectral(q1, params) + tauk;
    l1 = LL .* q1;

    % Third stage
    q2 = M .* (q + dt * (13861 * k0 / 62500 + 6889 * k1 / 62500 + 8611 * l0 / 62500 - 1743 * l1 / 31250));
    k2 = RHS_Spectral(q2, params) + tauk;
    l2 = LL .* q2;

    % Fourth stage
    q3 = M .* (q + dt * (-0.04884659515311858 * k0 - 0.1777206523264010 * k1 + 0.8465672474795196 * k2 + 0.1446368660269822 * l0 - 0.2239319076133447 * l1 + 0.4492950415863626 * l2));
    k3 = RHS_Spectral(q3, params) + tauk;
    l3 = LL .* q3;

    % Fifth stage
    q4 = M .* (q + dt * (-0.1554168584249155 * k0 - 0.3567050098221991 * k1 + 1.058725879868443 * k2 + 0.3033959883786719 * k3 + 0.09825878328356477 * l0 - 0.5915442428196704 * l1 + 0.8101210538282996 * l2 + 0.2831644057078060 * l3));
    k4 = RHS_Spectral(q4, params) + tauk;
    l4 = LL .* q4;

    % Sixth stage
    q5 = M .* (q + dt * (0.2014243506726763 * k0 + 0.008742057842904184 * k1 + 0.1599399570716811 * k2 + 0.4038290605220775 * k3 + 0.2260645738906608 * k4 + 0.1579162951616714 * l0 + 0.1867589405240008 * l2 + 0.6805652953093346 * l3 - 0.2752405309950067 * l4));
    k5 = RHS_Spectral(q5, params) + tauk;
    l5 = LL .* q5;

    % Successful step, proceed to evaluation
    q_next = real(ifft2(q + dt * (0.1579162951616714 * (k0 + l0) + 0.1867589405240008 * (k2 + l2) + 0.6805652953093346 * (k3 + l3) - 0.2752405309950067 * (k4 + l4) + (k5 + l5) / 4)));
    q_next = fft2(q_next);
end

function tauk = compute_tauk(params, q_size)
    % Compute the surface wind stress in spectral space
    [KX, KY] = meshgrid(0:q_size(1)-1, 0:q_size(2)-1);
    [X, Y] = meshgrid(-pi:2*pi/q_size(1):pi-2*pi/q_size(1), -pi:2*pi/q_size(2):pi-2*pi/q_size(2));

    if strcmp(params.forcingtype, 'sin')
        tau = params.tau0 * sin(Y).^2;
    elseif strcmp(params.forcingtype, 'sinh')
        tau = params.tau0 * sech(Y / params.tdel).^2;
    elseif strcmp(params.forcingtype, 'constant')
        tau = params.tau0;
    else
        tau = zeros(size(Y));
    end

    tauk = -1i * KY .* fft2(tau);
end

function name = casename(N, Ltracers, sigmax, Tend, spectrum_type, damping_type, dtsave)
    % Create a case name based on parameters
    name = sprintf('N%d_L%d_sigmax%.2f_Tend%.1f_spec%s_damp%s_dtsave%.2f', ...
                   N, Ltracers, sigmax, Tend, spectrum_type, damping_type, dtsave);
end
