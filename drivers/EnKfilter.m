%-------------------------------------------------------------------------%
% solver for lagrangian tracers and filter using barotropic data
%-------------------------------------------------------------------------%

% DO NOT SAVE COVARIANCE

    function [pin, m, Xt] = EnKfilter(pin)

        % Designing and implementing an ensemble Kalman filter (EnKF) for a Lagrangian data assimilation problem in MATLAB involves several steps. Here's a basic outline of the process:

        % Problem Formulation: 
        % Define the Lagrangian data assimilation problem, including the dynamical model, observational data, and the goal of assimilating observations into the model.
        
        % Ensemble Generation:
        % Generate an ensemble of model states that represent the uncertainty in the initial conditions. This ensemble will be used to estimate the state covariance.
        
        % Ensemble Prediction: 
        % Propagate each ensemble member forward in time using the dynamical model.
        
        % Observation Update:
        %  Update each ensemble member using observational data, adjusting the state estimates based on the observations.
        
        % Ensemble Covariance Update:
        %  Update the ensemble covariance matrix based on the observations and the ensemble mean.
        
        % Assimilation Loop:
        %  Repeat steps 3-5 for each assimilation cycle
        


% lets begin ::

%-------------------------------------------------------------------------%
% solver for lagrangian tracers and filter using barotropic data (modified for EnKF with recovered data)
%-------------------------------------------------------------------------%

% DO NOT SAVE COVARIANCE



    %-------------------------------------------------------------------------%
    % inputs
    %-------------------------------------------------------------------------%
    
    Tend = pin.Tend;
    
    damping_type = pin.damping_type;% d + mu |k|^2
    d = pin.d;
    if ~isfield(pin,'mu')
      pin.mu = 0;
    end
    mu = pin.mu;
    
    spectrum_type = pin.spectrum_type;
    if ~isfield(pin,'E0')
      pin.E0 = 0;
    end
    if ~isfield(pin,'num_modes_flat')
      pin.num_modes_flat = 0;
    end
    E0 = pin.E0;
    num_modes_flat = pin.num_modes_flat; % number of modes with flat spectrum
    
    
    % sigmax = pin.sigmax;
    sigmax = 0.1;
    Ltracers = pin.Ltracers;
    N = pin.N; % highest wavenumber
    dt = pin.dt;
    dtsave = pin.dtsave;
    
    % defaults
    
    
    if ~isfield(pin,'fsavedata')
      pin.fsavedata = false;
    end
    fsavedata = pin.fsavedata;
    
    
    if ~isfield(pin,'datafolder')
      pin.datafolder = '';
    end
    datafolder = pin.datafolder;
    %-------------------------------------------------------------------------%
    % setup
    %-------------------------------------------------------------------------%
    F = 0; % zero forcing or periodic forcing only (not tested)
    
    t = 0:dt:Tend+dt;
    n_t = length(t);
    
    dtsave = round(dtsave/dt)*dt; % round to a factor of dt
    nsavef = dtsave/dt;
    tsave = t(1:nsavef:length(t));
    n_tsave = length(tsave);
    
    
    fprintf('Total modes (2*N+1) = %i^2, max wavenumber resolved %d\n', N, N/2-1)
    Ntotal = N*N;
    
    L = 2*pi;
    x = L/N*(0:N-1)';  % x coordinate
    y = L/N*(0:N-1)';  % y coordinate
    [xx,yy] = meshgrid(x,y); % xy grid
    
    kx = [0:N/2-1 0 -N/2+1:-1]; % x wave vector
    ky = [0:N/2-1 0 -N/2+1:-1]; % y wave vector
    [kxx, kyy] = meshgrid(kx, ky);
    knorm = hypot(kxx,kyy);
    
    projx = (1i*-kyy./knorm);
    projy = (1i* kxx./knorm);
    projx(zeromodes2d(N)) = 0;
    projy(zeromodes2d(N)) = 0;
    
    if strcmp(damping_type, 'scaledep')
      gamma = (d + mu*knorm.^2);
    elseif strcmp(damping_type, 'unif')
      gamma = d*ones(N,N);
    elseif strcmp(damping_type,'exp')
      gamma = d;
    else
      error('unkown damping type')
    end
    
    if strcmp(spectrum_type, 'kolmo')
      sigma = sqrt(2*gamma.*(E0*knorm.^(-5/3))./knorm);
    elseif strcmp(spectrum_type, 'com_kolmo')
      sigma = sqrt(2*gamma.*(E0*(knorm/num_modes_flat).^(-5/3)));
      sigma(knorm < num_modes_flat) = sqrt(2*gamma(knorm < num_modes_flat)*E0./knorm(knorm < num_modes_flat));
    elseif strcmp(spectrum_type, 'third')
      sigma = sqrt(2*gamma.*(E0*knorm.^(-3))./knorm);
    elseif strcmp(spectrum_type, 'com_third')
        sigma = sqrt(2*gamma.*(E0*(knorm/num_modes_flat).^(-3)));
        sigma(knorm < num_modes_flat) = sqrt(2*gamma(knorm < num_modes_flat)*E0./knorm(knorm < num_modes_flat));
      elseif strcmp(spectrum_type, 'fourth')
        sigma = sqrt(2*gamma.*(E0*knorm.^(-4))./knorm);
      elseif strcmp(spectrum_type, 'flat')
        sigma = sqrt(2*gamma.*(E0*knorm.^(0))./knorm);
      elseif strcmp(spectrum_type, 'exp')
        if isscalar(pin.sigmax)
          sigma = pin.sigmax*ones(N,N);
        else
          sigma = pin.sigmax;
        end
      else
        error('unkown spectrum type')
      end
      
      sigma(junkmodes2d(N)) = 0;
      Sigma = sigma(:).*speye(Ntotal);
      gamma(junkmodes2d(N)) = 0;
      Gamma = gamma(:).' .* speye(Ntotal);
      
      
      %-------------------------------------------------------------------------%
      % save output parameter struct
      
      pout.tsave = tsave;
      pout.dtsave = dtsave;
      pout.N = N;
      pout.Ntotal = Ntotal;
      
      pout.sigma = sigma;
      pout.Sigma = Sigma;
      
      pout.gamma = gamma;
      pout.Gamma = Gamma;
      
      pout.x = x;
      pout.y = y;
      pout.kx = kx;
      pout.ky = ky;
      pout.kxx = kxx;
      pout.kyy = kyy;
      
      
      %-------------------------------------------------------------------------%
      % solve filter for velocity recovery ENKF
      %-------------------------------------------------------------------------%
      
      % R = zeros(N*N, N*N, n_tsave);
      m = zeros(N*N, n_tsave);
      X = pin.X; % Truth data (X) for EnKF
      
      
      % EnKF parameters (adjust as needed)
      ensemble_size = 2; % Number of ensemble members
      
      % Initialize ensemble with recovered data (m)
      m_ensemble = repmat(m, 1, ensemble_size); % Initial guess for ensemble
      is = 1;
      display(n_t); 
      for i = 1:n_t
    
        H = eye(N*N); 
        z = X(:, 1); 
        P_f = cov(m_ensemble'); 
        K = P_f * H' * inv(H * P_f * H' + sigmax^-2 * eye(N*N));
        
    for j = 1:ensemble_size
      if iscell(z)
          z_numeric = cell2mat(z);
      else
          z_numeric = z;
      end

      H_m_ensemble = H * m_ensemble(:, j);
      if iscell(H_m_ensemble)
          H_m_ensemble_numeric = cell2mat(H_m_ensemble);
      else
          H_m_ensemble_numeric = H_m_ensemble;
      end

      if ~isequal(size(z_numeric), size(H_m_ensemble_numeric))
          if numel(z_numeric) == numel(H_m_ensemble_numeric)
              if numel(z_numeric) == length(z_numeric)
                  z_numeric = reshape(z_numeric, [], 1);
                  H_m_ensemble_numeric = reshape(H_m_ensemble_numeric, [], 1);
              elseif numel(z_numeric) == length(H_m_ensemble_numeric)
                  z_numeric = reshape(z_numeric, 1, []);
                  H_m_ensemble_numeric = reshape(H_m_ensemble_numeric, 1, []);
              else
                  warning('Incompatible dimensions between z and H_m_ensemble');
                  continue; 
              end
          else
              warning('Incompatible dimensions between z and H_m_ensemble');
              continue;     
       end
        end

      % Update ensemble members
      m_ensemble(:, j) = m_ensemble(:, j) + K * (z_numeric - H_m_ensemble_numeric);
  end
       
    
        m_bar = mean(m_ensemble, 2);
        

        if mod(i, nsavef) == 0
          is = is + 1;
      %    R(:,:,is) = Rprev;
          m(:,is) = m_bar; 
        end
      end
    % fprintf('filter time %s\n', time2str(toc(timer)));
      
% Xt = X(:,1:nsavef:length(t)); 



pin.computed = pout; % save computed parameters

casestr = casename(pin.N, pin.Ltracers, sigmax, pin.Tend, pin.spectrum_type, pin.damping_type, pin.dtsave);
datafolder = fullfile(datafolder,casestr);
if fsavedata && ~isfolder(datafolder); mkdir(datafolder); end

if fsavedata
    save(fullfile(datafolder, 'data_enKF.mat'), 'm','pin', '-v7.3')
end



end 