function ensemble_updated = EnKF(ensemble, y, H, GAMMA, XI, model_func, dt)
    % Ensemble Kalman Filter Algorithm
    % input:
    %     ensemble: ensemble members at the previous time step
    %     y: current observation
    %     H: observation matrix (linear observation model)
    %     GAMMA: observation noise covariance matrix
    %     XI: dynamics noise covariance matrix
    %     model_func: dynamics model function
    %     dt: time step
    % output:
    %     ensemble_updated: updated ensemble members
    
    N_ENS = size(ensemble, 1);
    d_STATE = size(ensemble, 2);
    
    % Part 1: Prediction
    for i = 1:N_ENS
        ensemble(i, :) = model_func(ensemble(i, :), dt, XI);
    end
    % Compute predicted mean by averaging ensembles
    ens_mean = mean(ensemble, 1);
    ens_residual = ensemble - ens_mean;
    % Compute forecast error covariance matrix C
    C = cov(ens_residual);
    
    % Part 2: Analysis
    % Compute Kalman gain K
    S = H * C * H' + GAMMA;
    K = C * H' / S;
    ensemble_updated = zeros(size(ensemble));
    for i = 1:N_ENS
        % Sample eta from Normal(0, GAMMA)
        eta = mvnrnd(zeros(1, size(GAMMA, 1)), GAMMA);
        % Perturb the observation with the sampled eta
        y_mod = y + eta;
        % Compute the updated ensemble member
        v = (eye(d_STATE) - K * H) * ensemble(i, :)' + K * y_mod';
        ensemble_updated(i, :) = v';
    end
    
    return
end
