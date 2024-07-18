function [pin, m, Xt, ensemble] = lagrangianEnKF(pin)
    % Inputs and setup remain largely unchanged
    Tend = pin.Tend;
    damping_type = pin.damping_type;
    d = pin.d;
    mu = pin.mu;
    spectrum_type = pin.spectrum_type;
    E0 = pin.E0;
    num_modes_flat = pin.num_modes_flat;
    Ltracers = pin.Ltracers;
    N = pin.N;
    dt = pin.dt;
    dtsave = pin.dtsave;
    fsavedata = pin.fsavedata;
    datafolder = pin.datafolder;
    
    % Ensemble Kalman Filter specific setup
    Nens = 50; % Number of ensemble members
    ensemble = zeros(N*N, Nens); % Initialize ensemble matrix
    X = cell2mat(pin.X); % Assuming X contains your initial state
    for i = 1:Nens
        % Perturb initial conditions for each ensemble member
        ensemble(:,i) = X(:,1) + 0.1*randn(size(X,1),1); % Adjust perturbation as needed
    end
    
    % Your model dynamics and observation operator go here
    % model_dynamics = @(state, dt) ...;
    % observation_operator = @(state) ...;
    
    % Main simulation loop
    for i = 1:n_t
        % Forecast step: Propagate each ensemble member
        for j = 1:Nens
            ensemble(:,j) = model_dynamics(ensemble(:,j), dt);
        end
        
        % Update step: Only if observations are available
        if mod(i, observation_interval) == 0
            observations = get_observations_for_time(i*dt); % You need to implement this
            H = observation_operator; % Your observation operator
            R = cov(observations); % Observation error covariance
            for j = 1:Nens
                % Calculate Kalman Gain
                Hx = H * ensemble(:,j);
                P = cov(ensemble');
                K = P * H' / (H * P * H' + R);
                
                % Update ensemble member
                ensemble(:,j) = ensemble(:,j) + K * (observations - Hx);
            end
        end
        
        % Optional: Resampling or ensemble adjustment
        
        % Saving or processing results at intervals
        if mod(i, save_interval) == 0
            % Save or process as needed
        end
    end
    
    % Post-processing
    % Compute mean state, spread, etc. from ensemble
    % m = mean(ensemble, 2);
    % Xt = ...; % Extract or compute specific outputs as needed
    
    if fsavedata
        % Include saving mechanism for ensemble, mean state, etc.
        save(fullfile(datafolder, 'ensembleData.mat'), 'ensemble', 'm', 'Xt', '-v7.3');
    end
    
    end
    