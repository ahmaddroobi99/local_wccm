function [pin, m, Xt] = LDA_ENKF_QG_localization(pin)
    % Initialization
    [pin, params, x, y, kxx, kyy, projx, projy, tsave, nsavef, dtsave, N, Ntotal] = initialize_QG(pin);
    
    % EnKF with localization begins here
    ensemble_size = 35; % Increased ensemble size for better state representation
    mprev = zeros(N*N, 1);
    X_ens = mvnrnd(mprev', eye(N*N), ensemble_size)';
    m = zeros(N*N, numel(tsave)); % Array to store the results
    X = cell2mat(pin.X); % Observation matrix
    dXdt = diff(X, 1, 2) / pin.dt;
    dXdt = [dXdt, dXdt(:,end)]; % Extend dXdt to match the time dimension of X

    timer = tic; % Start timer
    Px = zeros(2 * pin.Ltracers, N * N); % Matrix for projections
    pin.nt= pin.Tend/pin.dt;
    
    % Main Loop
    for i = 1:pin.nt
        % Compute projection matrix Px
        Px = compute_projection_matrix(Px, kxx, kyy, projx, projy, X, pin.Ltracers, i);

        % EnKF update
        [X_ens, X_mean] = enkf_update(X_ens, dXdt(:, i), Px, pin.sigmax, ensemble_size, N);

        % QG Dynamics Forecast for the ensemble
        for j = 1:ensemble_size
            q_next = QG_dynamics(reshape(X_ens(:, j), N, N), pin.dt, params);
            X_ens(:, j) = reshape(q_next, [], 1);
        end

        X_ens_reshaped = mean(X_ens, 2);

        if mod(i, nsavef) == 0
            m(:, floor(i / nsavef) + 1) = X_ens_reshaped;
        end

        if any(isnan(m(:)))
            error('NaN in m detected at iteration %d!', i);
        end
    end

    fprintf('Filter time %s\n', time2str(toc(timer)));

    % Save only dt_save time steps
    Xt = X(:, 1:nsavef:length(pin.tlen));

    casestr = casename(pin.N, pin.Ltracers, pin.sigmax, pin.Tend, pin.spectrum_type, pin.damping_type, pin.dtsave);
    datafolder = fullfile(pin.datafolder, casestr);

    if pin.fsavedata && ~isfolder(datafolder)
        mkdir(datafolder);
    end

    if pin.fsavedata
        save(fullfile(datafolder, 'qgahmad_data.mat'), 'm', 'Xt', 'pin', '-v7.3');
    end
end

function [pin, params, x, y, kxx, kyy, projx, projy, tsave, nsavef, dtsave, N, Ntotal] = initialize_QG(pin)
    % Extract input parameters
    Tend = pin.Tend;
    damping_type = pin.damping_type;  
    d = pin.d;
    mu = pin.mu;
    spectrum_type = pin.spectrum_type;
    E0 = pin.E0;   
    num_modes_flat = 1;
    N = pin.N; 
    dt = pin.dt;
    dtsave = pin.dtsave;
    fsavedata = pin.fsavedata; 
    datafolder = pin.datafolder;
    t = 0:dt:Tend+dt;
    dtsave = round(dtsave/dt)*dt; % Round to a multiple of dt
    nsavef = dtsave/dt;
    tsave = t(1:nsavef:length(t));
    Ntotal = N*N;
    pin.tlen=length(t);



    % QG parameters
    qlim = 1.5E5; % if any q > qlim, simulation stops
    beta = 0.05;
    dd = 5E-3;
    epsf = 0;
    kf = 50; % no effect since epsf = 0
    sigf = 1; % no effect since epsf = 0
    alpf = 0; % no effect since epsf = 0
    tau0 = 0.04;
    tdel = 0; % no effect for sin forcing
    forcing_type = 'sin';
    Nt = 2e5;
    sh = 8; % hyperdiffusion exponent (needs to be even)
    nu = calculate_nu(N, kf, sh, epsf, tau0);

    

    % Setup grid and wave vectors
    L = 2*pi;
    x = L/N*(0:N-1)';
    y = L/N*(0:N-1)';
    [knorm,kx,ky,kxx, kyy, projx, projy] = setup_wave_vectors(N, x, y, d, mu,damping_type, spectrum_type, E0, num_modes_flat);

    gamma = calculate_gamma(N,d, mu, knorm, damping_type);
    sigma = calculate_sigma(N,gamma, knorm, E0, spectrum_type, num_modes_flat);
    sigma(junkmodes2d(N)) = 0;
    Sigma = sigma(:) .* speye(Ntotal);
    gamma(junkmodes2d(N)) = 0;
    Gamma = gamma(:) .* speye(Ntotal);


    params = struct('beta', beta, 'sh', sh, 'dd', dd, 'N', N, 'nu', nu, ...
                    'kf', kf, 'epsf', epsf, 'sigf', sigf, 'alpf', alpf, ...
                    'tau0', tau0, 'tdel', tdel, 'forcing_type', forcing_type);

    pin.computed = struct('tsave', tsave, 'dtsave', dtsave, 'N', N, 'Ntotal', Ntotal,'gamma',gamma,'sigma',sigma,'Sigma',Sigma,'Gamma',Gamma,'x', x, 'y', y, 'kx', kx, 'ky', ky, 'kxx', kxx, 'kyy', kyy);


end

function nu = calculate_nu(N, kf, sh, epsf, tau0)
    if N == 16
        nu = 2.5e-28 *(50/50)^(-2*sh+2/3);
    elseif N == 256
        nu = 5E-33 *(kf/50)^(-2*sh+2/3);
        if tau0 == 0
            nu = nu*100;
        end
    elseif N == 512
        nu = 1E-34 *(kf/50)^(-2*sh+2/3);
    end
end

function [knorm,kx,ky,kxx, kyy, projx, projy] = setup_wave_vectors(N, x, y, d, mu,damping_type, spectrum_type, E0, num_modes_flat)
    [xx, yy] = meshgrid(x, y);
    kx = [0:N/2-1 0 -N/2+1:-1];
    ky = [0:N/2-1 0 -N/2+1:-1];
    [kxx, kyy] = meshgrid(kx, ky);
    knorm = hypot(kxx, kyy);
    projx = (1i * -kyy ./ knorm);
    projy = (1i * kxx ./ knorm);
    projx(zeromodes2d(N)) = 0;
    projy(zeromodes2d(N)) = 0;
    
end

function gamma = calculate_gamma(N,d, mu, knorm, damping_type)
    if strcmp(damping_type, 'scaledep')
        gamma = (d + mu * knorm.^2);
    elseif strcmp(damping_type, 'unif')
        gamma = d * ones(N, N);
    else
        error('Unknown damping type')
    end
end

function sigma = calculate_sigma(N,gamma, knorm, E0, spectrum_type, num_modes_flat)
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
end

function Px = compute_projection_matrix(Px, kxx, kyy, projx, projy, X, Ltracers, i)
    for n = 1:Ltracers
        Xh = exp(1i*(kxx*X((n-1)*2 + 1, i) + kyy*X((n-1)*2 + 2, i)));
        x1 = Xh.*projx;
        x2 = Xh.*projy;
        Px((n-1)*2 + 1,:) = x1(:).';
        Px((n-1)*2 + 2,:) = x2(:).';
    end
end

function [X_ens, X_mean] = enkf_update(X_ens, dXdt_i, Px, sigmax, ensemble_size, N)
    X_mean = mean(X_ens, 2);  % Mean of the ensemble
    X_prime = X_ens - X_mean; % Deviations from the mean
    C = X_prime * X_prime' / (ensemble_size - 1); % Ensemble covariance

    % Apply localization
    localization_radius = 10; % Example localization radius
    localization_matrix = create_localization_matrix(N, localization_radius);
    C = C .* localization_matrix; % Apply localization

    % Add regularization to avoid singular matrix
    C = C + 1e-5 * eye(size(C));

    % Compute the Kalman gain
    Cov = Px * C * Px' + sigmax^2 * eye(size(Px, 1));  % Corrected observation noise variance
    K = C * Px' / Cov;

    % Apply the EnKF update
    X_ens_diff = dXdt_i - Px * X_ens;
    X_ens = X_ens + K * X_ens_diff;

    % Apply covariance inflation
    inflation_factor = 1.05; % Example inflation factor
    X_ens = X_mean + inflation_factor * (X_ens - X_mean);
end

function q_next = QG_dynamics(q, dt, params)
    N = params.N;
    tauk = compute_tauk(params, size(q));
    k = [0:N/2 -N/2+1:-1]'; % wavenumbers
    LL = compute_LL(N, k, params.nu, params.sh);
    M = 1 ./ (1 - .25 * dt * LL); % ARK4 stability function

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
end

function LL = compute_LL(N, k, nu, sh)
    LL = zeros(N, N);
    for jj = 1:N
        for ii = 1:N
            kr = sqrt(k(ii)^2 + k(jj)^2);
            LL(ii, jj) = -nu * kr^(2 * sh);
        end
    end
end

function tauk = compute_tauk(params, q_size)
    [KX, KY] = meshgrid(0:q_size(1)-1, 0:q_size(2)-1);
    [X, Y] = meshgrid(-pi:2*pi/q_size(1):pi-2*pi/q_size(1), -pi:2*pi/q_size(2):pi-2*pi/q_size(2));

    if strcmp(params.forcing_type, 'sin')
        tau = params.tau0 * sin(Y).^2;
    elseif strcmp(params.forcing_type, 'sinh')
        tau = params.tau0 * sech(Y / params.tdel).^2;
    elseif strcmp(params.forcing_type, 'constant')
        tau = params.tau0;
    else
        tau = zeros(size(Y));
    end

    tauk = -1i * KY .* fft2(tau);
end

function localization_matrix = create_localization_matrix(N, radius)
    localization_matrix = zeros(N*N, N*N);
    for i = 1:N
        for j = 1:N
            for k = 1:N
                for l = 1:N
                    distance = sqrt((i - k)^2 + (j - l)^2);
                    if distance <= radius
                        localization_matrix(sub2ind([N N], i, j), sub2ind([N N], k, l)) = 1;
                    end
                end
            end
        end
    end
end
