% this code testing nothing just see saff 




function [pin, X] = lagrangianfilter_baro(inputParams)

    % Function to implement the Lagrangian tracer filter using barotropic data
    %
    % INPUT PARAMETERS:
    %   inputParams          - A struct containing necessary input parameters
    %                           defined as follows:
    %       Tend              - Total simulation time
    %       dampingType       - Damping type used in the model (string)
    %       d                 - Value associated with the specified dampingType
    %       mu                - Additional damping term coefficient (optional, default=0)
    %       spectrumType      - Spectrum type used in the model (string)
    %       E0                - Energy scale for the energy spectrum
    %       numModesFlat      - Number of flat spectrum modes (optional, default=0)
    %       sigmax            - Standard deviation of measurement noise
    %       Ltracers          - Number of Lagrangian tracers
    %       N                 - Highest wavenumber resolution
    %       dt                - Time step size
    %       dtsave            - Sampling interval for saved data points
    %       fsavedata         - Flag indicating whether to store intermediate data
    %       dataFolder        - Path to folder where data files are stored (optional, empty string by default)
    
    % OUTPUT STRUCTURE:
    %   pin         - Struct containing necessary information about the solution
    %                           including:
    %       ts                - Array of simulated times
    %       tsSave            - Vector of sampled time indices for storage
    %       N                 - Size of the domain
    %       totalN            - Square of N
    %       sigma             - Matrix representing standard deviation of process noise
    %       gamma             - Matrix representing viscosity term
    %       x                 - Coordinate array along the x direction
    %       y                 - Coordinate array along the y direction
    %       kx                - Wavevector array along the x direction
    %       ky                - Wavevector array along the y direction
    %       knorm             - Magnitude of the 2D wavevectors
    %       projx             - Projection matrix onto x component
    %       projy             - Projection matrix onto y component
    %       U                 - Velocity field matrix
    %       V                 - Velocity field matrix
    %       phi               - Streamfunction matrix
    %       omega             - Vorticity matrix
    %       covarianceMatrix  - Covariance matrix of the state
    %       K                 - Kalman gain matrix
    %       Z                 - Measurement matrix
    %       X                 - State vector
    %       P                 - Error covariance matrix
    %       Q                 - Process noise covariance matrix
    %       R                 - Measurement noise covariance matrix
    %       KGainHistory      - History of Kalman gains over time
    %       estimatedVelocity - Estimated velocity matrix
    %       measuredVelocity  - Measured velocity matrix
    %       errVector         - Error between true and estimated velocities
    %       errMag            - Magnitude of the error vector
    %       innovHist         - Innovations history over time
    
    % Set default values
    if isnan(inputParams.mu) || isinf(inputParams.mu)
        inputParams.mu = 0;
    end
    if isnan(inputParams.numModesFlat) || isinf(inputParams.numModesFlat)
        inputParams.numModesFlat = 0;
    end
    
    % Define constants
    F = 0; % Zero forcing
    constPi = pi;
    
    % Calculate derived quantities
    t = 0:inputParams.dt:inputParams.Tend + inputParams.dt;
    nt = length(t);
    dtsaveRounded = round(inputParams.dtsave/inputParams.dt)*inputParams.dt;
    nSaveF = dtsaveRounded/inputParams.dt;
    tsave = t(1:nSaveF:length(t));
    nTSave = length(tsave);
    L = constPi;
    dx = L/inputParams.N;
    dy = dx;
    [XX YY] = ndgrid(dx*(0:inputParams.N-1)', dy*(0:inputParams.N-1)');
    KX = repmat(((0:inputParams.N/2-1)', zeros(1, inputParams.N/2))', inputParams.N, 1);
    KY = repmat(permute(KX, [2 1]), [1 inputParams.N/2]);
    KNorm = sqrt(KX.^2 + KY.^2) + eps;
    PROJX = -1i * KY ./ KNorm;
    PROJY = 1i * KX ./ KNorm;
    PROJX(KX == 0 & KY == 0) = 0;
    PROJY(KX == 0 & KY == 0) = 0;

    switch upper(inputParams.dampingType)
        case 'DAMPINGTYPE_SCALEDEP'
            GAMMA = inputParams.d + inputParams.mu .* KNorm.^2;
        case 'DAMPINGTYPE_UNIFORM'
            GAMMA = inputParams.d * ones(inputParams.N, inputParams.N);
        otherwise
            error('Invalid choice of damping type.');
    end

    switch upper(inputParams.spectrumType)
        case 'SPECTRUMTYPE_KOLMO'
            SIGMA = sqrt(2 * GAMMA .* (inputParams.E0 .* KNorm.^(-5/3)) ./ KNorm);
        case 'SPECTRUMTYPE_THIRD'
            SIGMA = sqrt(2 * GAMMA .* (inputParams.E0 .* KNorm.^(-3)) ./ KNorm);
        case 'SPECTRUMTYPE_FLAT'
            SIGMA = sqrt(2 * GAMMA .* (inputParams.E0 .* KNorm.^(0)) ./ KNorm);
        case 'SPECTRUMTYPE_EXP'
            SIGMA = inputParams.sigmax * ones(inputParams.N, inputParams.N);
        otherwise
            error('Invalid choice of spectrum type.');
    end

    SIGMA(KX == 0 & KY == 0) = 0;
    SIGMA_vec = diag(SIGMA(:));
    GAMMA_vec = diag(GAMMA(:));

    % Create structures for output
    pin.ts = t;
    pin.tsSave = tsave;
    pin.N = inputParams.N;
    pin.totalN = inputParams.N^2;
    pin.sigma = SIGMA;
    pin.gamma = GAMMA;
    pin.x = XX;
    pin.y = YY;
    pin.kx = KX;
    pin.ky = KY;
    pin.kNorm = KNorm;
    pin.projx = PROJX;
    pin.projy = PROJY;

      % Continue with the code...

    % Preallocate arrays
    U = zeros([size(XX) nTSave]);
    V = zeros([size(YY) nTSave]);
    PHI = zeros([size(YY) nTSave]);
    OMEGA = zeros([size(YY) nTSave]);
    COVARIANCE_MATRIX = zeros([inputParams.N^2 inputParams.N^2 nTSave]);
    K_GAIN_MATRIX = zeros([inputParams.N^2 inputParams.Ltracers nTSave]);
    Z = zeros([inputParams.Ltracers nTSave]);
    X = zeros([inputParams.N^2 1 nTSave]);
    P = zeros([inputParams.N^2 inputParams.N^2 nTSave]);
    Q = SIGMA_vec * SIGMA_vec;
    R = inputParams.sigmax^2 * eye(inputParams.Ltracers);
    K_GAIN_HISTORY = zeros([inputParams.N^2 inputParams.Ltracers nTSave]);
    ESTIMATED_VELOCITY = zeros([size(YY) nTSave]);
    MEASURED_VELOCITY = zeros([size(YY) nTSave]);
    ERR_VECTOR = zeros([size(YY) nTSave]);
    ERR_MAG = zeros([1 nTSave]);
    INNOV_HIST = zeros([inputParams.Ltracers nTSave]);

    % Initial conditions
    X(:, :, 1) = mvnrnd(zeros(inputParams.N^2, 1), SIGMA_vec);

    % Time-stepping loop
    for tIdx = 1:nTSave
        [u, v, phi, omega] = barotropic_model(inputParams.N, inputParams.d, inputParams.mu, F, GAMMA, SIGMA, PROJX, PROJY, X(:, :, tIdx), XX, YY, KX, KY);
        U(:, :, tIdx) = u;
        V(:, :, tIdx) = v;
        PHI(:, :, tIdx) = phi;
        OMEGA(:, :, tIdx) = omega;

        % Store intermediate data
        if inputParams.fsavedata
            COVARIANCE_MATRIX(:, :, tIdx) = reshape(cov(X(:, :, tIdx)), inputParams.N^2, inputParams.N^2);
            K_GAIN_MATRIX(:, :, tIdx) = X(:, :, tIdx) * (X(:, :, tIdx)' * X(:, :, tIdx) + R)^-1;
            Z(:, tIdx) = H * X(:, :, tIdx) + mvnrnd(zeros(inputParams.Ltracers, 1), R)';
            K_GAIN_HISTORY(:, :, tIdx) = K_GAIN_MATRIX(:, :, tIdx);
            ESTIMATED_VELOCITY(:, :, tIdx) = reshape(H * X(:, :, tIdx), inputParams.N, inputParams.N);
            MEASURED_VELOCITY(:, :, tIdx) = reshape(Z(:, tIdx), inputParams.N, inputParams.N);
            ERR_VECTOR(:, :, tIdx) = reshape(H * X(:, :, tIdx) - Z(:, tIdx), inputParams.N, inputParams.N);
            ERR_MAG(tIdx) = norm(ERR_VECTOR(:, :, tIdx), 'fro');
            INNOV_HIST(:, tIdx) = Z(:, tIdx) - H * X(:, :, tIdx);
        end

        % Time update (Prediction)
        X(:, :, tIdx + 1) = X(:, :, tIdx) + F * X(:, :, tIdx) + mvnrnd(zeros(inputParams.N^2, 1), Q)';
        P(:, :, tIdx + 1) = P(:, :, tIdx) + F * P(:, :, tIdx) * F' + Q;

        % Measurement update (Correction)
        K = P(:, :, tIdx + 1) * H' * inv(H * P(:, :, tIdx + 1) * H' + R);
        X(:, :, tIdx + 1) = X(:, :, tIdx + 1) + K * (Z(:, tIdx + 1) - H * X(:, :, tIdx + 1));
        P(:, :, tIdx + 1) = (eye(inputParams.N^2) - K * H) * P(:, :, tIdx + 1);

        % Damping update
        switch upper(inputParams.dampingType)
            case 'DAMPINGTYPE_SCALEDEP'
                GAMMA = inputParams.d + inputParams.mu * KNorm.^2;
            case 'DAMPINGTYPE_UNIFORM'
                GAMMA = inputParams.d * ones(inputParams.N, inputParams.N);
        end
    end

    % Assign output
    pin.U = U;
    pin.V = V;
    pin.phi = PHI;
    pin.omega = OMEGA;
    pin.covarianceMatrix = COVARIANCE_MATRIX;
    pin.kGainMatrix = K_GAIN_MATRIX;
    pin.Z = Z;
    pin.X = X;
    pin.P = P;
    pin.Q = Q;
    pin.R = R;
    pin.kGainHistory = K_GAIN_HISTORY;
    pin.estimatedVelocity = ESTIMATED_VELOCITY;
    pin.measuredVelocity = MEASURED_VELOCITY;
    pin.errVector = ERR_VECTOR;
    pin.errMag = ERR_MAG;
    pin.innovHist = INNOV_HIST;
end

function [u, v, phi, omega] = barotropic_model(N, d, mu, F, gamma, sigma, projx, projy, X, XX, YY, KX, KY)
    % Function to compute the barotropic model
    
    % Calculate velocity components in Fourier space
    Uhat = -1i * KX .* X;
    Vhat = -1i * KY .* X;
    
    % Inverse Fourier transform to get velocity components in physical space
    u = real(ifft2(Uhat));
    v = real(ifft2(Vhat));
    
    % Calculate streamfunction and vorticity
    phi = real(ifft2(-KX .* Vhat + KY .* Uhat));
    omega = real(ifft2(-KX .* Vhat - KY .* Uhat));
    
    % Apply damping
    u = u .* exp(-gamma);
    v = v .* exp(-gamma);
    phi = phi .* exp(-gamma);
    omega = omega .* exp(-gamma);
    
    % Add forcing
    u = u + F * XX;
    v = v + F * YY;
    phi = phi + F * XX;
    omega = omega + F * YY;
end
