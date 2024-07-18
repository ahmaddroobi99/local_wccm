function [pin, m, Xt] = Algorithm_design(pin)
    Ne = 100; 
    Np = 100; 

    state_ensemble = initialize_state_ensemble(Ne, pin);
    param_particles = initialize_param_particles(Np);

    save_steps = length(1:pin.dtsave:pin.Tend);
    m = zeros(size(state_ensemble, 2), save_steps);
    Xt = zeros(size(state_ensemble, 2), save_steps);

    % Main filtering loop
    save_idx = 1;
    for k = 1:length(pin.dt:pin.dt:pin.Tend)
        % State estimation using EnKF
        state_ensemble = propagate_state_ensemble(state_ensemble, pin.dt, pin);
        observations = get_observations(k, pin.X);
        state_ensemble = update_state_ensemble(state_ensemble, observations, pin);

        % Parameter estimation using PF
        param_particles = propagate_param_particles(param_particles, pin.dt);
        weights = evaluate_likelihood(param_particles, observations, state_ensemble, pin);
        param_particles = resample_particles(param_particles, weights);

        % Save results at specified intervals
        if mod(k, pin.dtsave / pin.dt) == 0
            m(:, save_idx) = mean(state_ensemble, 1); % Save mean state
            Xt(:, save_idx) = observations; % Save true state (or observations)
            save_idx = save_idx + 1;
        end
    end
end

% Helper Functions

function state_ensemble = initialize_state_ensemble(Ne, pin)
    state_ensemble = randn(Ne, pin.N); % Example initialization
end

function param_particles = initialize_param_particles(Np)
    param_particles.d = 0.1 + 0.01 * randn(Np, 1); % Example initialization
    param_particles.mu = 0.05 + 0.005 * randn(Np, 1); % Example initialization
end

function state_ensemble = propagate_state_ensemble(state_ensemble, dt, pin)
    for i = 1:size(state_ensemble, 1)
        state_ensemble(i, :) = model_dynamics(state_ensemble(i, :), dt, pin);
    end
end

function state_ensemble = update_state_ensemble(state_ensemble, observations, pin)
    % Ensure observations are numeric and a column vector
    if iscell(observations)
        observations = cell2mat(observations);
    end
    observations = observations(:); % Ensure it is a column vector

    % Calculate the mean and covariance of the ensemble
    x_mean = mean(state_ensemble, 1);
    P = cov(state_ensemble);

    % Update step
    H = eye(size(x_mean, 2)); % Observation matrix (identity for direct observations)
    R = pin.sigmax^2 * eye(size(H, 1)); % Observation noise covariance

    for i = 1:size(state_ensemble, 1)
        state_vector = state_ensemble(i, :)';
        % Ensure state_vector is numeric and a column vector
        if iscell(state_vector)
            state_vector = cell2mat(state_vector);
        end
        state_vector = state_vector(:); % Ensure it is a column vector

        K = P * H' / (H * P * H' + R); % Kalman gain
        state_ensemble(i, :) = (state_vector + K * (observations - H * state_vector))';
    end
end

function param_particles = propagate_param_particles(param_particles, dt)
    param_particles.d = param_particles.d + 0.001 * randn(size(param_particles.d));
    param_particles.mu = param_particles.mu + 0.0005 * randn(size(param_particles.mu));
end

function weights = evaluate_likelihood(param_particles, observations, state_ensemble, pin)
    % Ensure observations are numeric
    if iscell(observations)
        observations = cell2mat(observations);
    end

    % Calculate likelihoods
    likelihoods = zeros(size(param_particles.d));
    for i = 1:length(param_particles.d)
        % Example likelihood calculation (can be replaced with actual likelihood computation)
        predicted_observations = state_ensemble * param_particles.d(i); % Simplified model
        likelihoods(i) = exp(-0.5 * sum((observations - predicted_observations).^2) / pin.sigmax^2);
    end
    weights = likelihoods / sum(likelihoods);
end

function param_particles = resample_particles(param_particles, weights)
    indices = randsample(1:length(weights), length(weights), true, weights);
    param_particles.d = param_particles.d(indices);
    param_particles.mu = param_particles.mu(indices);
end

function observations = get_observations(k, X)
    observations = X(:, k); % Example observation extraction
    % Ensure observations are numeric
    if iscell(observations)
        observations = cell2mat(observations);
    end
end

function new_state = model_dynamics(current_state, dt, pin)
    % Linear Stochastic Model (LSM) dynamics

    % Define the linear deterministic part (e.g., advection matrix A)
    A = eye(pin.N); % This should be replaced with your actual model matrix A

    % Propagate the state using the deterministic part
    deterministic_part = current_state * A; 

    % Add stochastic noise
    process_noise = sqrt(dt) * randn(size(current_state)); % Adjust the noise level as needed
    new_state = deterministic_part + process_noise;
end
