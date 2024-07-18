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


% 
% 
% 
% function [pin, m, Xt] = L(pin)
% % LDA_ENKF_QG: Implementation of the Ensemble Kalman Filter (EnKF) for Lagrangian Data Assimilation in the Quasi-Geostrophic (QG) model.
% 
% % Extract input parameters
% Tend = pin.Tend;
% damping_type = pin.damping_type;  
% d = pin.d;
% mu = pin.mu;
% spectrum_type = pin.spectrum_type;
% E0 = pin.E0;
% num_modes_flat = 2;
% sigmax = 0.1; % Assumed constant value, can be parameterized as needed
% Ltracers = pin.Ltracers;
% N = pin.N; 
% dt = pin.dt;
% dtsave = pin.dtsave;
% fsavedata = pin.fsavedata; 
% datafolder = pin.datafolder;
% 
% % Setup simulation parameters
% t = 0:dt:Tend+dt;
% n_t = length(t);
% dtsave = round(dtsave/dt)*dt; % Round to a multiple of dt
% nsavef = dtsave/dt;
% tsave = t(1:nsavef:length(t));
% n_tsave = length(tsave);
% 
% fprintf('Total modes (2*N+1) = %i^2, max wavenumber resolved %d\n', N, N/2-1);
% Ntotal = N*N;
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% N = 8;
% beta = 0.05;
% dd = 5E-3;
% 
% epsf = 0;
% kf = 50;  % no effect since epsf = 0
% sigf = 1; % no effect since epsf = 0
% alpf = 0; % no effect since epsf = 0
% tau0 = 0.04;
% tdel = 0; % no effect for sin forcing
% forcing_type = 'sin';
% 
% Nt = 2e5;
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % QG parameters
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% qlim = 1.5E5; % if any q > qlim, simulation stops
% 
% % Set physical parameters
% % epsf forcing amplitude
% % dd  % linear Ekman damping
% 
% % tau0  % surface wind stress strength
% % tdel  % width of wind stress
% 
% sh = 8; % hyperdiffusion exponent (needs to be even)
% if N == 8
%     if epsf == 0
%         nu = 2.5e-28 *(50/50)^(-2*sh+2/3); %5E-33; % Coefficient of biharmonic vorticity diffusion
%     else
%         nu = 2.5e-28 *(kf/50)^(-2*sh+2/3); %5E-33; % Coefficient of biharmonic vorticity diffusion
%     end
% elseif N == 256
%     nu = 5E-33 *(kf/50)^(-2*sh+2/3); %5E-33; % Coefficient of biharmonic vorticity diffusion
%     if tau0 == 0
%         nu = nu*100;
%     end
% elseif N == 512
%     nu = 1E-34 *(kf/50)^(-2*sh+2/3);
% end
% 
% 
% % Put useful stuff into a struct
% params = struct('beta',beta, 'sh',sh, 'dd',dd, 'nu',nu, 'N',N,...
%     'kf',kf, 'epsf',epsf, 'sigf',sigf, 'alpf',alpf,...
%     'tau0',tau0, 'tdel',tdel,'forcingtype',forcing_type);
% 
% 
% 
% 
% F = 0; % zero forcing or periodic forcing only (not tested)
% 
% t = 0:dt:Tend+dt;
% n_t = length(t);
% 
% dtsave = round(dtsave/dt)*dt; % round to a factor of dt
% nsavef = dtsave/dt;
% tsave = t(1:nsavef:length(t));
% n_tsave = length(tsave);
% 
% 
% fprintf('Total modes (2*N+1) = %i^2, max wavenumber resolved %d\n', N, N/2-1)
% Ntotal = N*N;
% 
% L = 2*pi;
% x = L/N*(0:N-1)';             % x coordinate
% y = L/N*(0:N-1)';             % y coordinate
% [xx,yy] = meshgrid(x,y);      % xy grid
% 
% kx = [0:N/2-1 0 -N/2+1:-1]; % x wave vector
% ky = [0:N/2-1 0 -N/2+1:-1]; % y wave vector
% [kxx, kyy] = meshgrid(kx, ky);
% knorm = hypot(kxx,kyy);
% 
% projx = (1i*-kyy./knorm);
% projy = (1i* kxx./knorm);
% projx(zeromodes2d(N)) = 0;
% projy(zeromodes2d(N)) = 0;
% 
% if strcmp(damping_type, 'scaledep')
%     gamma = (d + mu*knorm.^2);
% elseif strcmp(damping_type, 'unif')
%     gamma = d*ones(N,N);
% elseif strcmp(damping_type,'exp')
%     gamma = d;
% else
%     error('unkown damping type')
% end
% 
% if strcmp(spectrum_type, 'kolmo')
%     sigma = sqrt(2*gamma.*(E0*knorm.^(-5/3))./knorm);
% elseif strcmp(spectrum_type, 'com_kolmo')
%     sigma = sqrt(2*gamma.*(E0*(knorm/num_modes_flat).^(-5/3)));
%     sigma(knorm < num_modes_flat) = sqrt(2*gamma(knorm < num_modes_flat)*E0./knorm(knorm < num_modes_flat));
% elseif strcmp(spectrum_type, 'third')
%     sigma = sqrt(2*gamma.*(E0*knorm.^(-3))./knorm);
% elseif strcmp(spectrum_type, 'com_third')
%     sigma = sqrt(2*gamma.*(E0*(knorm/num_modes_flat).^(-3)));
%     sigma(knorm < num_modes_flat) = sqrt(2*gamma(knorm < num_modes_flat)*E0./knorm(knorm < num_modes_flat));
% elseif strcmp(spectrum_type, 'fourth')
%     sigma = sqrt(2*gamma.*(E0*knorm.^(-4))./knorm);
% elseif strcmp(spectrum_type, 'flat')
%     sigma = sqrt(2*gamma.*(E0*knorm.^(0))./knorm);
% elseif strcmp(spectrum_type, 'exp')
%     if isscalar(pin.sigma)
%         sigma = pin.sigma*ones(N,N);
%     else
%         sigma = pin.sigma;
%     end
% else
%     error('unkown spectrum type')
% end
% 
% sigma(junkmodes2d(N)) = 0;
% Sigma = sigma(:).*speye(Ntotal);
% gamma(junkmodes2d(N)) = 0;
% Gamma = gamma(:).*speye(Ntotal);
% 
% 
% % Save computed parameters
% pout = struct('tsave', tsave, 'dtsave', dtsave, 'N', N, 'Ntotal', Ntotal, ...
%     'sigma', sigma, 'Sigma', Sigma, 'gamma', gamma, 'Gamma', Gamma, ...
%     'x', x, 'y', y, 'kx', kx, 'ky', ky, 'kxx', kxx, 'kyy', kyy);
% pin.computed = pout;
% 
% % Initial ensemble setup
% % Initial ensemble setup
% ensemble_size = 10;
% mprev = zeros(N*N, 1);
% X_ens = mvnrnd(mprev', eye(N*N), ensemble_size)';
% m = zeros(N*N, n_tsave); % Array to store the results
% X = cell2mat(pin.X);     % Observation matrix
% dXdt = [(X(:,2)-X(:,1)) (X(:,2:end) - X(:,1:end-1))] / dt; % Time derivative of observations
% 
% timer = tic;             % Start timer
% Px = zeros(2 * Ltracers, N * N); % Matrix for projections
% % Sigma2 = Sigma * Sigma'; % Covariance of the model error
% 
% % Initial ensemble setup
% q = zeros(N, N, ensemble_size); % Initialize q for the ensemble with zeros
% 
% % Main loop
% for i = 1:n_t
%     % Compute projection matrix Px
%     for n = 1:Ltracers
%         Xh = exp(1i * (kxx * X((n-1)*2 + 1, i) + kyy * X((n-1)*2 + 2, i)));
%         Px((n-1)*2 + 1, :) = reshape(Xh .* projx, 1, []);
%         Px((n-1)*2 + 2, :) = reshape(Xh .* projy, 1, []);
%     end
% 
%     % Ensemble Kalman Filter (EnKF) update
%     X_mean = mean(X_ens, 2);  % Mean of the ensemble
%     X_prime = X_ens - X_mean; % Deviations from the mean
%     C = X_prime * X_prime' / (ensemble_size - 1); % Ensemble covariance
% 
%     % Compute the Kalman gain
%     Cov = Px * C * Px' + sigmax^-2 * eye(size(Px, 1)); 
%     K = C * Px' / Cov;
% 
%     % Surface wind stress
%     k = [0:N/2 -N/2+1:-1]'; % Wavenumbers
%     dx = 2 * pi / N;
%     [KX, KY] = meshgrid(k, k);
%     [X, Y] = meshgrid(-pi:dx:pi-dx, -pi:dx:pi-dx);
% 
%     if strcmp(forcing_type, 'sin')
%         tau = tau0 * sin(Y).^2;
%     elseif strcmp(forcing_type, 'sinh')
%         tau = tau0 * sech(Y/tdel).^2;
%     elseif strcmp(forcing_type, 'constant')
%         tau = tau0;
%     end
%     tauk = -1i * KY .* fft2(tau);
% 
%     tol = 1E-1;
%     r0 = 0.8 * tol;
% 
%     % QG Dynamics Forecast for the ensemble
%     for j = 1:ensemble_size
%         % Perform QG dynamics update for q using the QG solver
% 
%         L = zeros(N,N);
% 
%         M = 1./(1 - .25 * dt * L);
% 
%         % First stage ARK4
%         k0 = RHS_Spectral(q(:,:,j), params) + tauk;
%         l0 = L .* q(:,:,j);
% 
%         % Second stage
% 
%         q1 = M .* (q(:,:,j) + .5 * dt * k0 + .25 * dt * l0);
%         k1 = RHS_Spectral(q1, params) + tauk;
%         l1 = L .* q1;
% 
%         % Third stage
% 
%         q2 = M .* (q(:,:,j) + dt * (13861 * k0 / 62500 + 6889 * k1 / 62500 + 8611 * l0 / 62500 - 1743 * l1 / 31250));
%         k2 = RHS_Spectral(q2, params) + tauk;
%         l2 = L .* q2;
% 
%         % Fourth stage
% 
%         q3 = M .* (q(:,:,j) + dt * (-0.04884659515311858 * k0 - 0.1777206523264010 * k1 + 0.8465672474795196 * k2 + 0.1446368660269822 * l0 - 0.2239319076133447 * l1 + 0.4492950415863626 * l2));
%         k3 = RHS_Spectral(q3, params) + tauk;
%         l3 = L .* q3;
% 
%         % Fifth stage
% 
%         q4 = M .* (q(:,:,j) + dt * (-0.1554168584249155 * k0 - 0.3567050098221991 * k1 + 1.058725879868443 * k2 + 0.3033959883786719 * k3 + 0.09825878328356477 * l0 - 0.5915442428196704 * l1 + 0.8101210538282996 * l2 + 0.2831644057078060 * l3));
%         k4 = RHS_Spectral(q4, params) + tauk;
%         l4 = L .* q4;
% 
%         % Sixth stage
% 
%         q5 = M .* (q(:,:,j) + dt * (0.2014243506726763 * k0 + 0.008742057842904184 * k1 + 0.1599399570716811 * k2 + 0.4038290605220775 * k3 + 0.2260645738906608 * k4 + 0.1579162951616714 * l0 + 0.1867589405240008 * l2 + 0.6805652953093346 * l3 - 0.2752405309950067 * l4));
%         k5 = RHS_Spectral(q5, params) + tauk;
%         l5 = L .* q5;
% 
%         % Error control
% 
%         r1 = dt * max(max(max(abs(ifft2(0.003204494398459 * (k0 + l0) - 0.002446251136679 * (k2 + l2) - 0.021480075919587 * (k3 + l3) + 0.043946868068572 * (k4 + l4) - 0.023225035410765 * (k5 + l5))))));
% 
%         if r1 > tol
%             dt = .75 * dt;
%             continue;
%         end
% 
%         % Successful step, proceed to evaluation
%         t = t + dt;
%         qp = real(ifft2(q(:,:,j) + dt * (0.1579162951616714 * (k0 + l0) + 0.1867589405240008 * (k2 + l2) + 0.6805652953093346 * (k3 + l3) - 0.2752405309950067 * (k4 + l4) + (k5 + l5) / 4)));
% 
%         q(:,:,j) = fft2(qp);
%         dt = ((.75 * tol / r1)^(.3 / 4)) * ((r0 / r1)^(.4 / 4)) * dt;
%         r0 = r1;
%     end
% 
%     % Reshape q back to vector form
%     X_ens(:, j) = reshape(q(:,:,j), [], 1);
% 
%     % Apply the EnKF update
%     X_ens_diff = dXdt(:, i) - Px * X_ens;
%     X_ens_update = X_ens + K * X_ens_diff;
% 
%     % Update the EnKF ensemble
%     X_ens = X_ens_update;
%     X_ens_reshaped = mean(X_ens, 2);
% 
%     if mod(i, nsavef) == 0
%         m(:, floor(i / nsavef) + 1) = X_ens_reshaped;
%     end
% end
% fprintf('Filter time %s\n', time2str(toc(timer)))
% 
% % Save only dt_save time steps
% Xt = X(:, 1:nsavef:length(t)); 
% 
% pin.computed = pout; % Save computed parameters
% 
% casestr = casename(pin.N, pin.Ltracers, sigmax, pin.Tend, pin.spectrum_type, pin.damping_type, pin.dtsave);
% datafolder = fullfile(datafolder, casestr);
% 
% if fsavedata && ~isfolder(datafolder)
%     mkdir(datafolder);
% end
% 
% if fsavedata
%     save(fullfile(datafolder, 'hello_data.mat'), 'm', 'Xt', 'pin', '-v7.3')
% end
% end
% 
% 
% 
% 
% % % % % % % % % % % % % % % % % % % % % ----------------------------------------------
% function gamma = get_gamma(damping_type, d, mu, knorm)
%     if strcmp(damping_type, 'scaledep')
%         gamma = d + mu * knorm.^2;
%     elseif strcmp(damping_type, 'unif')
%         gamma = d * ones(size(knorm));
%     else
%         error('Unknown damping type');
%     end
% end
% 
% function sigma = get_sigma(spectrum_type, E0, knorm, gamma, num_modes_flat, N, sigmax)
%     if strcmp(spectrum_type, 'kolmo')
%         sigma = sqrt(2 * gamma .* (E0 * knorm.^(-5/3)) ./ knorm);
%     elseif strcmp(spectrum_type, 'com_kolmo')
%         sigma = sqrt(2 * gamma .* (E0 * (knorm / num_modes_flat).^(-5/3)) ./ knorm);
%     elseif strcmp(spectrum_type, 'unif')
%         sigma = sigmax * ones(size(knorm));
%     elseif strcmp(spectrum_type, 'none')
%         sigma = zeros(size(knorm));
%     else
%         error('Unknown spectrum type');
% % _________________________
% 
% 
% 
%     end
% end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% %  help :-> build an ROM fromQG itself like LSM
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% 
% 
% function [pin, m, Xt] = LDA_ENKF_QG(pin)
%     % LDA_ENKF_QG: Implementation of the Ensemble Kalman Filter (EnKF) for Lagrangian Data Assimilation in the Quasi-Geostrophic (QG) model.
% 
%     % Extract input parameters
%     Tend = pin.Tend;
%     damping_type = pin.damping_type;  
%     d = pin.d;
%     mu = pin.mu;
%     spectrum_type = pin.spectrum_type;
%     E0 = pin.E0;
%     num_modes_flat = 2;
%     sigmax = 0.1; % Assumed constant value, can be parameterized as needed
%     Ltracers = pin.Ltracers;
%     N = pin.N; 
%     dt = pin.dt;
%     dtsave = pin.dtsave;
%     fsavedata = pin.fsavedata; 
%     datafolder = pin.datafolder;
% 
%     % Setup simulation parameters
%     t = 0:dt:Tend+dt;
%     n_t = length(t);
%     dtsave = round(dtsave/dt)*dt; % Round to a multiple of dt
%     nsavef = dtsave/dt;
%     tsave = t(1:nsavef:length(t));
%     n_tsave = length(tsave);
% 
%     fprintf('Total modes (2*N+1) = %i^2, max wavenumber resolved %d\n', N, N/2-1);
%     Ntotal = N*N;
% 
%     % QG parameters
%     qlim = 1.5E5; % if any q > qlim, simulation stops
%     N = 8;
%     beta = 0.05;
%     dd = 5E-3;
% 
%     epsf = 0;
%     kf = 50;  % no effect since epsf = 0
%     sigf = 1; % no effect since epsf = 0
%     alpf = 0; % no effect since epsf = 0
%     tau0 = 0.04;
%     tdel = 0; % no effect for sin forcing
%     forcing_type = 'sin';
% 
%     Nt = 2e5;
% 
%     % Set physical parameters
%     sh = 8; % hyperdiffusion exponent (needs to be even)
%     if N == 8
%         if epsf == 0
%             nu = 2.5e-28 *(50/50)^(-2*sh+2/3); %5E-33; % Coefficient of biharmonic vorticity diffusion
%         else
%             nu = 2.5e-28 *(kf/50)^(-2*sh+2/3); %5E-33; % Coefficient of biharmonic vorticity diffusion
%         end
%     elseif N == 256
%         nu = 5E-33 *(kf/50)^(-2*sh+2/3); %5E-33; % Coefficient of biharmonic vorticity diffusion
%         if tau0 == 0
%             nu = nu*100;
%         end
%     elseif N == 512
%         nu = 1E-34 *(kf/50)^(-2*sh+2/3);
%     end
% 
% 
% 
%     % Put useful stuff into a struct
%     params = struct('beta', beta, 'sh', sh, 'dd', dd, 'nu', nu, 'N', N, ...
%                     'kf', kf, 'epsf', epsf, 'sigf', sigf, 'alpf', alpf, ...
%                     'tau0', tau0, 'tdel', tdel, 'forcingtype', forcing_type);
% 
%     F = 0; % zero forcing or periodic forcing only (not tested)
% 
%     % Setup grid and wave vectors
%     L = 2*pi;
%     x = L/N*(0:N-1)'; % x coordinate
%     y = L/N*(0:N-1)'; % y coordinate
%     [xx, yy] = meshgrid(x, y); % xy grid
% 
%     kx = [0:N/2-1 0 -N/2+1:-1]; % x wave vector
%     ky = [0:N/2-1 0 -N/2+1:-1]; % y wave vector
%     [kxx, kyy] = meshgrid(kx, ky);
%     knorm = hypot(kxx, kyy);
% 
%     projx = (1i * -kyy ./ knorm);
%     projy = (1i * kxx ./ knorm);
%     projx(zeromodes2d(N)) = 0;
%     projy(zeromodes2d(N)) = 0;
% 
%     % Determine gamma based on damping type
%     if strcmp(damping_type, 'scaledep')
%         gamma = (d + mu * knorm.^2);
%     elseif strcmp(damping_type, 'unif')
%         gamma = d * ones(N, N);
%     elseif strcmp(damping_type, 'exp')
%         gamma = d;
%     else
%         error('Unknown damping type')
%     end
% 
%     % Determine sigma based on spectrum type
%     if strcmp(spectrum_type, 'kolmo')
%         sigma = sqrt(2 * gamma .* (E0 * knorm.^(-5/3)) ./ knorm);
%     elseif strcmp(spectrum_type, 'com_kolmo')
%         sigma = sqrt(2 * gamma .* (E0 * (knorm / num_modes_flat).^(-5/3)));
%         sigma(knorm < num_modes_flat) = sqrt(2 * gamma(knorm < num_modes_flat) * E0 ./ knorm(knorm < num_modes_flat));
%     elseif strcmp(spectrum_type, 'third')
%         sigma = sqrt(2 * gamma .* (E0 * knorm.^(-3)) ./ knorm);
%     elseif strcmp(spectrum_type, 'com_third')
%         sigma = sqrt(2 * gamma .* (E0 * (knorm / num_modes_flat).^(-3)));
%         sigma(knorm < num_modes_flat) = sqrt(2 * gamma(knorm < num_modes_flat) * E0 ./ knorm(knorm < num_modes_flat));
%     elseif strcmp(spectrum_type, 'fourth')
%         sigma = sqrt(2 * gamma .* (E0 * knorm.^(-4)) ./ knorm);
%     elseif strcmp(spectrum_type, 'flat')
%         sigma = sqrt(2 * gamma .* (E0 * knorm.^(0)) ./ knorm);
%     elseif strcmp(spectrum_type, 'exp')
%         if isscalar(pin.sigma)
%             sigma = pin.sigma * ones(N, N);
%         else
%             sigma = pin.sigma;
%         end
%     else
%         error('Unknown spectrum type')
%     end
% 
%     sigma(junkmodes2d(N)) = 0;
%     Sigma = sigma(:) .* speye(Ntotal);
%     gamma(junkmodes2d(N)) = 0;
%     Gamma = gamma(:) .* speye(Ntotal);
% 
%     % Save computed parameters
%     pout = struct('tsave', tsave, 'dtsave', dtsave, 'N', N, 'Ntotal', Ntotal, ...
%                   'sigma', sigma, 'Sigma', Sigma, 'gamma', gamma, 'Gamma', Gamma, ...
%                   'x', x, 'y', y, 'kx', kx, 'ky', ky, 'kxx', kxx, 'kyy', kyy);
%     pin.computed = pout;
% 
%     % Initial ensemble setup
%     ensemble_size = 100;
%     mprev = zeros(N*N, 1);
%     X_ens = mvnrnd(mprev', eye(N*N), ensemble_size)';
%     m = zeros(N*N, n_tsave); % Array to store the results
%     X = cell2mat(pin.X);     % Observation matrix
%     dXdt = [(X(:,2)-X(:,1)) (X(:,2:end) - X(:,1:end-1))] / dt; % Time derivative of observations
% 
%     timer = tic;             % Start timer
%     Px = zeros(2 * Ltracers, N * N); % Matrix for projections
%     Sigma2 = Sigma * Sigma'; % Covariance of the model error
% 
%     % Main Loop
%     for i = 1:n_t
%         % Compute projection matrix Px
%         for n = 1:Ltracers
%             Xh = exp(1i * (kxx * X((n-1)*2 + 1, i) + kyy * X((n-1)*2 + 2, i)));
%             Px((n-1)*2 + 1, :) = reshape(Xh .* projx, 1, []);
%             Px((n-1)*2 + 2, :) = reshape(Xh .* projy, 1, []);
%         end
% 
%         % Ensemble Kalman Filter (EnKF) update
%         X_mean = mean(X_ens, 2);  % Mean of the ensemble
%         X_prime = X_ens - X_mean; % Deviations from the mean
%         C = X_prime * X_prime' / (ensemble_size - 1); % Ensemble covariance
% 
%         % Compute the Kalman gain
%         Cov = Px * C * Px' + sigmax^-2 * eye(size(Px, 1)); 
%         K = C * Px' / Cov;
% 
%         % Initialize q
%         if i == 1
%             q = fft2(zeros(N, N)); % Initial potential vorticity field in Fourier space
%         end
% 
%         % Update q using QG dynamics (example implementation)
%         for j = 1:ensemble_size
%             % Reshape q back to vector form for the current ensemble member
%             q_ens = reshape(X_ens(:, j), N, N);
% 
%             % Perform QG dynamics update for q 
%             % using the QG solver
%             % Example update step 
%             % (user should replace with 
%             % actual QG dynamics)
%             q_ens = QG_Step(q_ens, params, dt); % User-defined QG dynamics function
% 
%             % Reshape q back to vector form
%             X_ens(:, j) = q_ens(:);
%         end
% 
%         % Apply the EnKF update
%         X_ens_diff = dXdt(:, i) - Px * X_ens;
%         X_ens_update = X_ens + K * X_ens_diff;
% 
%         % Update the EnKF ensemble
%         X_ens = X_ens_update;
%         X_ens_reshaped = mean(X_ens, 2);
% 
%         if mod(i, nsavef) == 0
%             m(:, floor(i/nsavef) + 1) = X_ens_reshaped;
%         end
%     end
% 
%     fprintf('Filter time %s\n', time2str(toc(timer)))
% 
%     % Save only dt_save time steps
%     Xt = X(:, 1:nsavef:length(t)); 
% 
%     pin.computed = pout; % Save computed parameters
% 
%     casestr = casename(pin.N, pin.Ltracers, sigmax, pin.Tend, pin.spectrum_type, pin.damping_type, pin.dtsave);
%     datafolder = fullfile(datafolder, casestr);
% 
%     if fsavedata && ~isfolder(datafolder)
%         mkdir(datafolder);
%     end
% 
%     if fsavedata
%         save(fullfile(datafolder, 'hello_data.mat'), 'm', 'Xt', 'pin', '-v7.3')
%     end
% end
% 
% function q = QG_Step(q, params, dt)
%     % Example function to update the potential 
%     % vorticity field q using QG dynamics
% 
% 
%     % This is an ROM  and should 
%     % be replaced with actual QG dynamics
% 
%     % User-defined QG dynamics 
%     % update (ROM)
%     q = q + dt * (-params.beta * q + params.nu * laplacian(q));
% end
% 
% function L = laplacian(q)
%     % Example function to compute the Laplacian of q
%     % This is a placeholder and should be replaced with actual Laplacian computation
%     L = del2(q);
% end
% 
% 
% function str = time2str(time)
%     % Function to convert time to string (placeholder)
%     str = sprintf('%0.2f seconds', time);
% end


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%  now back to full QG dynamics   (Version 3)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


% function [pin, m, Xt] = LDA_ENKF_QG(pin)
%     % LDA_ENKF_QG: Implementation of the Ensemble Kalman Filter (EnKF) for Lagrangian Data Assimilation in the Quasi-Geostrophic (QG) model.
% 
%     % Extract input parameters
%     Tend = pin.Tend;
%     damping_type = pin.damping_type;  
%     d = pin.d;
%     mu = pin.mu;
%     spectrum_type = pin.spectrum_type;
%     E0 = pin.E0;
%     num_modes_flat = 2;
%     sigmax = 0.1; % Assumed constant value, can be parameterized as needed
%     Ltracers = pin.Ltracers;
%     N = pin.N; 
%     dt = pin.dt;
%     dtsave = pin.dtsave;
%     fsavedata = pin.fsavedata; 
%     datafolder = pin.datafolder;
% 
%     % Setup simulation parameters
%     t = 0:dt:Tend+dt;
%     n_t = length(t);
%     dtsave = round(dtsave/dt)*dt; % Round to a multiple of dt
%     nsavef = dtsave/dt;
%     tsave = t(1:nsavef:length(t));
%     n_tsave = length(tsave);
% 
%     fprintf('Total modes (2*N+1) = %i^2, max wavenumber resolved %d\n', N, N/2-1);
%     Ntotal = N*N;
% 
%     % QG parameters
%     qlim = 1.5E5; % if any q > qlim, simulation stops
%     beta = 0.05;
%     dd = 5E-3;
% 
%     epsf = 0;
%     kf = 50;  % no effect since epsf = 0
%     sigf = 1; % no effect since epsf = 0
%     alpf = 0; % no effect since epsf = 0
%     tau0 = 0.04;
%     tdel = 0; % no effect for sin forcing
%     forcing_type = 'sin';
% 
%     Nt = 2e5;
%     % Set physical parameters
%     sh = 8; % hyperdiffusion exponent (needs to be even)
%     if N == 8
%         if epsf == 0
%             nu = 2.5e-28 *(50/50)^(-2*sh+2/3); %5E-33; % Coefficient of biharmonic vorticity diffusion
%         else
%             nu = 2.5e-28 *(kf/50)^(-2*sh+2/3); %5E-33; % Coefficient of biharmonic vorticity diffusion
%         end
%     elseif N == 256
%         nu = 5E-33 *(kf/50)^(-2*sh+2/3); %5E-33; % Coefficient of biharmonic vorticity diffusion
%         if tau0 == 0
%             nu = nu*100;
%         end
%     elseif N == 512
%         nu = 1E-34 *(kf/50)^(-2*sh+2/3);
%     end
% 
%     % Put useful stuff into a struct
%     params = struct('beta', beta, 'sh', sh, 'dd', dd, 'nu', nu, 'N', N, ...
%                     'kf', kf, 'epsf', epsf, 'sigf', sigf, 'alpf', alpf, ...
%                     'tau0', tau0, 'tdel', tdel, 'forcingtype', forcing_type);
% 
%     F = 0; % zero forcing or periodic forcing only (not tested)
% 
%     % Setup grid and wave vectors
%     L = 2*pi;
%     x = L/N*(0:N-1)'; % x coordinate
%     y = L/N*(0:N-1)'; % y coordinate
%     [xx, yy] = meshgrid(x, y); % xy grid
% 
%     kx = [0:N/2-1 0 -N/2+1:-1]; % x wave vector
%     ky = [0:N/2-1 0 -N/2+1:-1]; % y wave vector
%     [kxx, kyy] = meshgrid(kx, ky);
%     knorm = hypot(kxx, kyy);
% 
%     projx = (1i * -kyy ./ knorm);
%     projy = (1i * kxx ./ knorm);
%     projx(zeromodes2d(N)) = 0;
%     projy(zeromodes2d(N)) = 0;
% 
%     % Determine gamma based on damping type
%     if strcmp(damping_type, 'scaledep')
%         gamma = (d + mu * knorm.^2);
%     elseif strcmp(damping_type, 'unif')
%         gamma = d * ones(N, N);
%     elseif strcmp(damping_type, 'exp')
%         gamma = d;
%     else
%         error('Unknown damping type')
%     end
% 
%     % Determine sigma based on spectrum type
%     if strcmp(spectrum_type, 'kolmo')
%         sigma = sqrt(2 * gamma .* (E0 * knorm.^(-5/3)) ./ knorm);
%     elseif strcmp(spectrum_type, 'com_kolmo')
%         sigma = sqrt(2 * gamma .* (E0 * (knorm / num_modes_flat).^(-5/3)));
%         sigma(knorm < num_modes_flat) = sqrt(2 * gamma(knorm < num_modes_flat) * E0 ./ knorm(knorm < num_modes_flat));
%     elseif strcmp(spectrum_type, 'third')
%         sigma = sqrt(2 * gamma .* (E0 * knorm.^(-3)) ./ knorm);
%     elseif strcmp(spectrum_type, 'com_third')
%         sigma = sqrt(2 * gamma .* (E0 * (knorm / num_modes_flat).^(-3)));
%         sigma(knorm < num_modes_flat) = sqrt(2 * gamma(knorm < num_modes_flat) * E0 ./ knorm(knorm < num_modes_flat));
%     elseif strcmp(spectrum_type, 'fourth')
%         sigma = sqrt(2 * gamma .* (E0 * knorm.^(-4)) ./ knorm);
%     elseif strcmp(spectrum_type, 'flat')
%         sigma = sqrt(2 * gamma .* (E0 * knorm.^(0)) ./ knorm);
%     elseif strcmp(spectrum_type, 'exp')
%         if isscalar(pin.sigma)
%             sigma = pin.sigma * ones(N, N);
%         else
%             sigma = pin.sigma;
%         end
%     else
%         error('Unknown spectrum type')
%     end
% 
%     sigma(junkmodes2d(N)) = 0;
%     Sigma = sigma(:) .* speye(Ntotal);
%     gamma(junkmodes2d(N)) = 0;
%     Gamma = gamma(:) .* speye(Ntotal);
% 
%     % Save computed parameters
%     pout = struct('tsave', tsave, 'dtsave', dtsave, 'N', N, 'Ntotal', Ntotal, ...
%                   'sigma', sigma, 'Sigma', Sigma, 'gamma', gamma, 'Gamma', Gamma, ...
%                   'x', x, 'y', y, 'kx', kx, 'ky', ky, 'kxx', kxx, 'kyy', kyy);
%     pin.computed = pout;
% 
%     % Initial ensemble setup
%     ensemble_size = 10;
%     mprev = zeros(N*N, 1);
%     X_ens = mvnrnd(mprev', eye(N*N), ensemble_size)';
%     m = zeros(N*N, n_tsave); % Array to store the results
%     X = cell2mat(pin.X);     % Observation matrix
%     dXdt = [(X(:,2)-X(:,1)) (X(:,2:end) - X(:,1:end-1))] / dt; % Time derivative of observations
% 
%     timer = tic;             % Start timer
%     Px = zeros(2 * Ltracers, N * N); % Matrix for projections
% 
%     % Main Loop
%     for i = 1:n_t
%         % Compute projection matrix Px
%         for n = 1:Ltracers
%             Xh = exp(1i * (kxx * X((n-1)*2 + 1, i) + kyy * X((n-1)*2 + 2, i)));
%             Px((n-1)*2 + 1, :) = reshape(Xh .* projx, 1, []);
%             Px((n-1)*2 + 2, :) = reshape(Xh .* projy, 1, []);
%         end
% 
%         % Ensemble Kalman Filter (EnKF) update
%         X_mean = mean(X_ens, 2);  % Mean of the ensemble
%         X_prime = X_ens - X_mean; % Deviations from the mean
%         C = X_prime * X_prime' / (ensemble_size - 1); % Ensemble covariance
% 
%         % Compute the Kalman gain
%         Cov = Px * C * Px' + sigmax^-2 * eye(size(Px, 1)); 
%         K = C * Px' / Cov;
% 
%         % Initialize q
%         if i == 1
%             q = fft2(zeros(N, N)); % Initial potential vorticity field in Fourier space
%         end
% 
%         % Update q using QG dynamics
%         for j = 1:ensemble_size
%             % Reshape q back to vector form for the current ensemble member
%             q_ens = reshape(X_ens(:, j), N, N);
% 
%             % Perform QG dynamics update for q using the QG solver
%             dtt = 1E-5; % initial time step size
%             q_ens = QG_Step(q_ens, params, dtt); % User-defined QG dynamics function
% 
%             % Reshape q back to vector form
%             X_ens(:, j) = q_ens(:);
%         end
% 
%         % Apply the EnKF update
%         X_ens_diff = dXdt(:, i) - Px * X_ens;
%         X_ens_update = X_ens + K * X_ens_diff;
% 
%         % Update the EnKF ensemble
%         X_ens = X_ens_update;
%         X_ens_reshaped = mean(X_ens, 2);
% 
%         if mod(i, nsavef) == 0
%             m(:, floor(i/nsavef) + 1) = X_ens_reshaped;
%         end
%     end
% 
%     fprintf('Filter time %s\n', time2str(toc(timer)))
% 
%     % Save only dt_save time steps
%     Xt = X(:, 1:nsavef:length(t)); 
% 
%     pin.computed = pout; % Save computed parameters
% 
%     casestr = casename(pin.N, pin.Ltracers, sigmax, pin.Tend, pin.spectrum_type, pin.damping_type, pin.dtsave);
%     datafolder = fullfile(datafolder, casestr);
% 
%     if fsavedata && ~isfolder(datafolder)
%         mkdir(datafolder);
%     end
% 
%     if fsavedata
%         save(fullfile(datafolder, 'hello_data.mat'), 'm', 'Xt', 'pin', '-v7.3')
%     end
% end
% 
% function q = QG_Step(q, params, dt)
%     % QG_Step: Update the potential vorticity field q using QG dynamics
%     % This function performs a single time step using ARK4 method
%     N=8;
%     L = zeros(N,N);
% 
%     % Coefficients for the ARK4 method
%     a21 = 0.5; a31 = 13861/62500; a32 = 6889/62500;
%     a41 = -0.04884659515311858; a42 = -0.1777206523264010; a43 = 0.8465672474795196;
%     a51 = -0.1554168584249155; a52 = -0.3567050098221991; a53 = 1.058725879868443; a54 = 0.3033959883786719;
%     a61 = 0.2014243506726763; a62 = 0.008742057842904184; a63 = 0.1599399570716811; a64 = 0.4038290605220775; a65 = 0.2260645738906608;
% 
%     b1 = 0.1579162951616714; b3 = 0.1867589405240008; b4 = 0.6805652953093346; b5 = -0.2752405309950067; b6 = 0.25;
% 
%     % Initialize variables
%     k0 = zeros(size(q));
%     k1 = zeros(size(q));
%     k2 = zeros(size(q));
%     k3 = zeros(size(q));
%     k4 = zeros(size(q));
%     k5 = zeros(size(q));
%     l0 = zeros(size(q));
%     l1 = zeros(size(q));
%     l2 = zeros(size(q));
%     l3 = zeros(size(q));
%     l4 = zeros(size(q));
%     l5 = zeros(size(q));
% 
%     % Compute the surface wind stress in spectral space
%     tauk = compute_tauk(params, size(q));
% 
%     % First stage ARK4
%     M = 1./(1 - 0.25 * dt * L);
%     k0 = RHS_Spectral(q, params) + tauk;
%     l0 = L .* q;
% 
%     % Second stage
%     q1 = M .* (q + a21 * dt * k0 + 0.25 * dt * l0);
%     k1 = RHS_Spectral(q1, params) + tauk;
%     l1 = L .* q1;
% 
%     % Third stage
%     q2 = M .* (q + dt * (a31 * k0 + a32 * k1) + dt * (8611/62500 * l0 - 1743/31250 * l1));
%     k2 = RHS_Spectral(q2, params) + tauk;
%     l2 = L .* q2;
% 
%     % Fourth stage
%     q3 = M .* (q + dt * (a41 * k0 + a42 * k1 + a43 * k2) + dt * (1446368660269822 * l0 - 2239319076133447 * l1 + 4492950415863626 * l2));
%     k3 = RHS_Spectral(q3, params) + tauk;
%     l3 = L .* q3;
% 
%     % Fifth stage
%     q4 = M .* (q + dt * (a51 * k0 + a52 * k1 + a53 * k2 + a54 * k3) + dt * (9825878328356477 * l0 - 5915442428196704 * l1 + 8101210538282996 * l2 + 2831644057078060 * l3));
%     k4 = RHS_Spectral(q4, params) + tauk;
%     l4 = L .* q4;
% 
%     % Sixth stage
%     q5 = M .* (q + dt * (a61 * k0 + a62 * k1 + a63 * k2 + a64 * k3 + a65 * k4) + dt * (1579162951616714 * l0 + 1867589405240008 * l2 + 6805652953093346 * l3 - 2752405309950067 * l4));
%     k5 = RHS_Spectral(q5, params) + tauk;
%     l5 = L .* q5;
% 
%     % Error control
%     tol = 1E-1;
%     r0 = 0.8 * tol;
%     r1 = dt * max(max(max(abs(ifft2(0.003204494398459 * (k0 + l0) - 0.002446251136679 * (k2 + l2))))));
%     if r1 > tol, dt = 0.75 * dt; return, end
% 
%     % Successful step, proceed to evaluation
%     q = real(ifft2(q + dt * (b1 * (k0 + l0) + b3 * (k2 + l2) + b4 * (k3 + l3) - b5 * (k4 + l4) + b6 * (k5 + l5))));
% 
%     % Update time step size for next iteration
%     dt = ((0.75 * tol / r1)^(0.3 / 4)) * ((r0 / r1)^(0.4 / 4)) * dt;
%     r0 = r1;
% end
% 
% function RHS = RHS_Spectral(q_hat, p)
%     % Function takes Fourier coefficients of PV (q_hat) and struct containing
%     % parameters (p) and evaluates RHS of barotropic QG equations except for
%     % high-k dissipation. Returns Fourier coefficients of RHS.
%     % Jacobian is not dealiased because solutions are well resolved.
% 
%     persistent dX dY DX DY Laplacian InvBT
%     if isempty(DX)
%         k = [0:p.N/2 -p.N/2+1:-1]';
%         dX = 1i * repmat(k', [p.N 1]);
%         dY = 1i * repmat(k, [1 p.N]);
%         Laplacian = dX.^2 + dY.^2;
%         InvBT = 1 ./ Laplacian; InvBT(1,1) = 0;
% 
%         k = [0:p.N/2-1 0 -p.N/2+1:-1]';
%         dX = 1i * repmat(k', [p.N 1]);
%         dY = 1i * repmat(k, [1 p.N]);
%         % For the dealiased jacobian:
%         k = [0:.75*p.N-1 0 -.75*p.N+1:-1]';
%         DX = 1i * repmat(k', [1.5*p.N 1]);
%         DY = 1i * repmat(k, [1 1.5*p.N]);
%         clear k
%     end
% 
%     % Invert for psi
%     q_bt = q_hat; %+hk(:,:);
%     psi_hat = InvBT .* q_hat;
% 
%     % Calculate Ekman plus beta plus mean shear
%     RHS = zeros([p.N p.N]);
%     RHS = -p.beta * dX .* psi_hat - p.dd * q_bt;
% 
%     % For using a 3/2-rule dealiased jacobian:
%     % physical space, 3/2 grid; factor of (9/4) scales fft
%     Psi_hat = zeros([1.5*p.N 1.5*p.N]);
%     Psi_hat(1:p.N/2+1,1:p.N/2+1) = (9/4) * psi_hat(1:p.N/2+1,1:p.N/2+1);
%     Psi_hat(1:p.N/2+1,p.N+2:1.5*p.N) = (9/4) * psi_hat(1:p.N/2+1,p.N/2+2:p.N);
%     Psi_hat(p.N+2:1.5*p.N,1:p.N/2+1) = (9/4) * psi_hat(p.N/2+2:p.N,1:p.N/2+1);
%     Psi_hat(p.N+2:1.5*p.N,p.N+2:1.5*p.N) = (9/4) * psi_hat(p.N/2+2:p.N,p.N/2+2:p.N);
%     Q_hat = zeros([1.5*p.N 1.5*p.N]);
%     Q_hat(1:p.N/2+1,1:p.N/2+1) = (9/4) * q_bt(1:p.N/2+1,1:p.N/2+1);
%     Q_hat(1:p.N/2+1,p.N+2:1.5*p.N) = (9/4) * q_bt(1:p.N/2+1,p.N/2+2:p.N);
%     Q_hat(p.N+2:1.5*p.N,1:p.N/2+1) = (9/4) * q_bt(p.N/2+2:p.N,1:p.N/2+1);
%     Q_hat(p.N+2:1.5*p.N,p.N+2:1.5*p.N) = (9/4) * q_bt(p.N/2+2:p.N,p.N/2+2:p.N);
%     % calculate u.gradq on 3/2 grid
%     u = real(ifft2(-DY .* Psi_hat));
%     v = real(ifft2(DX .* Psi_hat));
%     qx = real(ifft2(DX .* Q_hat));
%     qy = real(ifft2(DY .* Q_hat));
%     jaco_real = u .* qx + v .* qy;
%     % fft, 3/2 grid; factor of (4/9) scales fft
%     Jaco_hat = (4/9) * fft2(jaco_real);
%     % reduce to normal grid
%     jaco_hat = zeros([p.N p.N]);
%     jaco_hat(1:p.N/2+1,1:p.N/2+1) = Jaco_hat(1:p.N/2+1,1:p.N/2+1);
%     jaco_hat(1:p.N/2+1,p.N/2+2:p.N) = Jaco_hat(1:p.N/2+1,p.N+2:1.5*p.N);
%     jaco_hat(p.N/2+2:p.N,1:p.N/2+1) = Jaco_hat(p.N+2:1.5*p.N,1:p.N/2+1);
%     jaco_hat(p.N/2+2:p.N,p.N/2+2:p.N) = Jaco_hat(p.N+2:1.5*p.N,p.N+2:1.5*p.N);
% 
%     % Put it all together
%     RHS = RHS - jaco_hat;
% end
% 
% function tauk = compute_tauk(params, q_size)
%     % Compute the surface wind stress in spectral space
%     [KX, KY] = meshgrid(0:q_size(1)-1, 0:q_size(2)-1);
%     [X, Y] = meshgrid(-pi:2*pi/q_size(1):pi-2*pi/q_size(1), -pi:2*pi/q_size(2):pi-2*pi/q_size(2));
% 
%     if strcmp(params.forcingtype, 'sin')
%         tau = params.tau0 * sin(Y).^2;
%     elseif strcmp(params.forcingtype, 'sinh')
%         tau = params.tau0 * sech(Y / params.tdel).^2;
%     elseif strcmp(params.forcingtype, 'constant')
%         tau = params.tau0;
%     else
%         tau = zeros(size(Y));
%     end
% 
%     tauk = -1i * KY .* fft2(tau);
% end




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
    sigmax = 0.1; % Assumed constant value, can be  parameterized as needed
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
    params = struct('beta', beta, 'sh', sh, 'dd', dd, 'nu', nu, 'N', N, ...
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
    elseif strcmp(damping_type, 'exp')
        gamma = d;
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

    % Initial ensemble setup
    ensemble_size = 10;
    mprev = zeros(N*N, 1);
    X_ens = mvnrnd(mprev', eye(N*N), ensemble_size)';
    m = zeros(N*N, n_tsave); % Array to store the results
    X = cell2mat(pin.X);     % Observation matrix
    dXdt = [(X(:,2)-X(:,1)) (X(:,2:end) - X(:,1:end-1))]./dt;

    timer = tic;             % Start timer
    Px = zeros(2 * Ltracers, N * N); % Matrix for projections
    dt = 1E-5; % initial time step size


    % Initialize
    t = 0;
    qp(:,:) = 0.01*randn(params.N);
    qp(:,:) = qp(:,:) - mean(mean(qp(:,:)));
    q = fft2(qp);

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
            % Perform QG dynamics update for q using the QG solver

            % L = compute_L(N); % Define your operator L here


            
            L = zeros(N,N);
            tauk = compute_tauk(params,  size(q));

            M = 1./(1 - .25 * dt * L);

            % First stage ARK4
            k0 = RHS_Spectral(q, params) + tauk;
            l0 = L .* q;

            % Second stage
            q1 = M .* (q + .5 * dt * k0 + .25 * dt * l0);
            k1 = RHS_Spectral(q1, params) + tauk;
            l1 = L .* q1;

            % Third stage
            q2 = M .* (q + dt * (13861 * k0 / 62500 + 6889 * k1 / 62500 + 8611 * l0 / 62500 - 1743 * l1 / 31250));
            k2 = RHS_Spectral(q2, params) + tauk;
            l2 = L .* q2;

            % Fourth stage
            q3 = M .* (q + dt * (-0.04884659515311858 * k0 - 0.1777206523264010 * k1 + 0.8465672474795196 * k2 + 0.1446368660269822 * l0 - 0.2239319076133447 * l1 + 0.4492950415863626 * l2));
            k3 = RHS_Spectral(q3, params) + tauk;
            l3 = L .* q3;

            % Fifth stage
            q4 = M .* (q + dt * (-0.1554168584249155 * k0 - 0.3567050098221991 * k1 + 1.058725879868443 * k2 + 0.3033959883786719 * k3 + 0.09825878328356477 * l0 - 0.5915442428196704 * l1 + 0.8101210538282996 * l2 + 0.2831644057078060 * l3));
            k4 = RHS_Spectral(q4, params) + tauk;
            l4 = L .* q4;

            % Sixth stage
            q5 = M .* (q + dt * (0.2014243506726763 * k0 + 0.008742057842904184 * k1 + 0.1599399570716811 * k2 + 0.4038290605220775 * k3 + 0.2260645738906608 * k4 + 0.1579162951616714 * l0 + 0.1867589405240008 * l2 + 0.6805652953093346 * l3 - 0.2752405309950067 * l4));
            k5 = RHS_Spectral(q5, params) + tauk;
            l5 = L .* q5;

            % Error control
            r1 = dt * max(max(max(abs(ifft2(0.003204494398459 * (k0 + l0) - 0.002446251136679 * (k2 + l2) - 0.021480075919587 * (k3 + l3) + 0.043946868068572 * (k4 + l4) - 0.023225035410765 * (k5 + l5))))));
            tol = 1E-1;

            if r1 > tol
                dt = .75 * dt;
                continue;
            end

            % Successful step, proceed to evaluation
            t = t + dt;
            qp = real(ifft2(q + dt * (0.1579162951616714 * (k0 + l0) + 0.1867589405240008 * (k2 + l2) + 0.6805652953093346 * (k3 + l3) - 0.2752405309950067 * (k4 + l4) + (k5 + l5) / 4)));

            q = fft2(qp);
            % r0 = 0.8 * tol;
            % 
            % dt = ((.75 * tol / r1)^(.3 / 4)) * ((r0 / r1)^(.4 / 4)) * dt;
            % r0 = r1;

            % Reshape q back to vector form
            X_ens(:, j) = reshape(q, [], 1);
        end

           if any(isnan(X_ens(:)))
          fprintf('NaN is m Recovery matrix  solution detected!\n')
           end

        % Apply the EnKF update
        X_ens_diff = dXdt(:, i) - Px * X_ens;
        X_ens_update = X_ens + K * X_ens_diff;

        % Update the EnKF ensemble
        X_ens = X_ens_update;
        X_ens_reshaped = mean(X_ens, 2);

        if mod(i, nsavef) == 0
            m(:, floor(i / nsavef) + 1) = X_ens_reshaped;
        end

        if any(isnan(m(:)))
          fprintf('NaN is m Recovery matrix  solution detected!\n')
        end
    end

    fprintf('Filter time %s\n', time2str(toc(timer)));

    % Save only dt_save time steps
    Xt = X(:, 1:nsavef:length(t)); 

    pin.computed = pout; % Save computed parameters

    casestr = casename(pin.N, pin.Ltracers, sigmax, pin.Tend, pin.spectrum_type, pin.damping_type, pin.dtsave);
    datafolder = fullfile(datafolder, casestr);

    if fsavedata && ~isfolder(datafolder)
        mkdir(datafolder);
    end

    if fsavedata
        save(fullfile(datafolder, 'qqqqqqq_data.mat'), 'm', 'Xt', 'pin', '-v7.3')
    end
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
