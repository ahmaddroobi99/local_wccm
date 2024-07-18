%-------------------------------------------------------------------------%
% Ensemble Kalman Filter for velocity recovery
%-------------------------------------------------------------------------%

% Parameters
N = 100;  % State dimension
Ltracers = 10;  % Number of tracers
n_ens = 100;  % Ensemble size
n_t = 100;  % Number of time steps
nsavef = 10;  % Save frequency
dt = 0.1;  % Time step
sigmax = 1;  % Measurement noise
Sigma = diag(sigmax^2*ones(N,1));  % Measurement noise covariance
Gamma = diag(0.1*ones(N,1));  % Process noise covariance

% Initialize ensemble
X_ens = cell2mat(pin.X);  % Initialize ensemble
n_ens = 10;  % Ensemble size
X_ens = repmat(X_ens, 1, n_ens);  % Repeat for each ensemble member

dXdt = [(X_ens(:,2)-X_ens(:,1)) (X_ens(:,2:end) - X_ens(:,1:end-1))]./dt;

% initial conditions
is = 1;
Rprev = zeros(N*N, N*N);  %1/2*(Sigma*Sigma')./Gamma;
mprev = zeros(N*N,1);
% R(:,:,is) = Rprev;
m(:,is) = mprev;

timer = tic;
Px = zeros(2*Ltracers, N*N);

Sigma2 = Sigma*Sigma';

for i = 1:n_t
    for n = 1:Ltracers
        Xh = exp(1i*(kxx*X_ens((n-1)*2 + 1,i) + kyy*X_ens((n-1)*2 + 2,i)));
        x1 = Xh.*projx;
        x2 = Xh.*projy;
        Px((n-1)*2 + 1,:) = x1(:).';
        Px((n-1)*2 + 2,:) = x2(:).';
    end
    
    % Forecast step
    X_ens = forecast_model(X_ens);
    
    % Compute ensemble mean and covariance
    X_mean = mean(X_ens, 2);
    P_forecast = cov(X_ens');
    
    % Update step
    for j = 1:n_ens
        % Update equation
        X_ens(:,j) = update_equation(X_ens(:,j), measurement);
    end
    
    % Compute analysis ensemble
    X_analysis = X_ens;
    
    % Save the analysis ensemble
    if mod(i, nsavef) == 0
        is = is + 1;
        X_analysis_save(:,:,is) = X_analysis;
    end
end

fprintf('filter time %s\n', time2str(toc(timer)));


%-------------------------------------------------------------------------%
% save only dt_save time steps

Xt = X_analysis(:,1:nsavef:length(t)); 


% Forecast function
function X_forecast = forecast_model(X_current)
    A = diag(0.9*ones(N,1));  % System matrix
    X_forecast = A * X_current;
end

% Update function
function X_update = update_equation(X_forecast, measurement)
    H = ones(1,N);  % Measurement matrix
    K = X_forecast * H' / (H * X_forecast * H' + sigmax^2);
    X_update = X_forecast + K * (measurement - H * X_forecast);
end

