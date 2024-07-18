function X = generate_eulerian_obs(datafilename, Ltracers, sigmax, varargin)

    p = inputParser;
    addParameter(p, 'savedata', false, @islogical)
    
    parse(p, varargin{:});
    fsavedata = p.Results.savedata;
    
    
    fnameqh = fullfile(datafilename, 'qh11.bin');
    fnamet = fullfile(datafilename, 't11.bin');
    dataext = load(fullfile(datafilename, 'diag11.mat'));
    
    fidqh = fopen(fnameqh,'r');
    fidt = fopen(fnamet,'r');
    
    t = fread(fidt,'double');
    dt = t(2) - t(1);
    nt = length(t);
    N =  dataext.params.N;
    Ntotal = N*N;
    
    
    k = [0:N/2 -N/2+1:-1]';
    [KX, KY] = meshgrid(k,k);
    Knorm = hypot(KX,KY);
    Kinvnorm = 1./Knorm; Kinvnorm(1,1) = 0;
    dX = 1i*repmat(k', [N 1]);
    dY = 1i*repmat(k, [1 N]);
    Laplacian = dX.^2 + dY.^2;
    InvBT = 1./Laplacian; InvBT(1,1) = 0;
    
    L = 2*pi;
    x = L/N*(0:N-1)';             % x coordinate
    y = L/N*(0:N-1)';             % y coordinate
    [xx,yy] = meshgrid(x,y);      % xy grid
    
    % hi 2024/3/22 
    X = zeros(2*Ltracers, nt);  % Initialize Eulerian observations
    % X(:,:,1) = abs(ifft2(fft2(randn(N,N)) .* fft2(sigmax*randn(N,N))));
    X(:,1) = x(end)*rand(2*Ltracers,1);
    X = mat2cell(X, 2*ones(Ltracers,1), nt);
    timer = tic;
    for i = 1:(nt-1)
        if ~mod(i,100)
           percentdone = 100 * i / (nt-1);
           msg = sprintf('Percent done: %3.1f\n', percentdone);
           fprintf(msg);
        end
        
        qq = fread(fidqh,2*Ntotal,'double');
        qq = complex(qq(1:Ntotal), qq(Ntotal+1:2*Ntotal));
    
        qh = reshape(qq,N,N);
        psi_hat = InvBT.*qh;
    
        v = real(ifft2(dX.*psi_hat));
        u = real(ifft2(-dY.*psi_hat));
        
        % % Update Eulerian observations
        % X(:,:,i+1) = abs(ifft2(fft2(X(:,:,i)) + fft2(u.*dt) + fft2(sigmax*randn(N,N)*sqrt(dt))));
        % X(:,:,i+1) = min(max(X(:,:,i+1),0),1); % Ensure observations are within valid range
        


        % so lets work in Eulerian

        % parfor n = 1:Ltracers
        %     Xu = mod(X{n}(1,i), x(end));
        %     Xv = mod(X{n}(2,i), x(end));
            
        %     uxi = interp2(xx, yy, u, Xu, Xv);
        %     vxi = interp2(xx, yy, v, Xu, Xv);
            
        %     X{n}(1,i+1) = X{n}(1,i) + uxi*dt + sigmax*randn()*sqrt(dt);
        %     X{n}(2,i+1) = X{n}(2,i) + vxi*dt + sigmax*randn()*sqrt(dt);
        % end

            % so lets work in Eulerian

    for n = 1:Ltracers
        Xu = mod(X{n}(1,1), x(end)); % Use the initial x-coordinate instead of the current one
        Xv = mod(X{n}(2,1), y(end)); % Use the initial y-coordinate instead of the current one
        
        uxi = interp2(xx, yy, u, Xu, Xv);
        vxi = interp2(xx, yy, v, Xu, Xv);
        
        X{n}(:,i+1) = X{n}(:,i); % Keep the spatial coordinates fixed
        X{n}(1,i+1) = X{n}(1,1) + uxi*dt + sigmax*randn()*sqrt(dt); % Update only time index
        X{n}(2,i+1) = X{n}(2,1) + vxi*dt + sigmax*randn()*sqrt(dt); % Update only time index

    end




    end
    fprintf('\n tracer solve time %s\n', time2str(toc(timer)));
    
    if fsavedata
        casestr = sprintf('Eul_obs4_Ltracers%i_sigmax%.2g', Ltracers, sigmax);
        save(fullfile(datafilename, strcat(casestr,'.mat')), 'X', 'sigmax', 'Ltracers', '-v7.3')
    end
    
end

% function X = generate_baro_obs(datafilename, Ltracers, sigmax, varargin)

%     p = inputParser;
%     addParameter(p, 'savedata', false, @islogical)
    
%     parse(p, varargin{:});
%     fsavedata = p.Results.savedata;
    
    
%     fnameqh = fullfile(datafilename, 'qh11.bin');
%     fnamet = fullfile(datafilename, 't11.bin');
%     dataext = load(fullfile(datafilename, 'diag11.mat'));
    
%     fidqh = fopen(fnameqh,'r');
%     fidt = fopen(fnamet,'r');
    
%     t = fread(fidt,'double');
%     dt = t(2) - t(1);
%     nt = length(t);
%     N =  dataext.params.N;
%     Ntotal = N*N;
    
    
%     k = [0:N/2 -N/2+1:-1]';
%     [KX, KY] = meshgrid(k,k);
%     Knorm = hypot(KX,KY);
%     Kinvnorm = 1./Knorm; Kinvnorm(1,1) = 0;
%     dX = 1i*repmat(k', [N 1]);
%     dY = 1i*repmat(k, [1 N]);
%     Laplacian = dX.^2 + dY.^2;
%     InvBT = 1./Laplacian; InvBT(1,1) = 0;
    
%     L = 2*pi;
%     x = L/N*(0:N-1)';             % x coordinate
%     y = L/N*(0:N-1)';             % y coordinate
%     [xx, yy] = meshgrid(x, y);    % xy grid
    
%     % grid_spacing = floor(sqrt(2 * Ltracers)); % Adjust this value as needed
%     % [X_grid, Y_grid] = meshgrid(1:grid_spacing:N, 1:grid_spacing:N);
    
%     % X_init = [X_grid(:), Y_grid(:)];
    
%     % if size(X_init, 1) < 2 * Ltracers
%     %     X_init = repmat(X_init, ceil(2 * Ltracers / size(X_init, 1)), 1);
%     % end
    
%     % X = mat2cell(X_init(1:(2 * Ltracers), :)', 2, ones(1, 2 * Ltracers));



%     X = zeros(2*Ltracers, nt);
%     X(:,1) = x(end)*rand(2*Ltracers,1);
%     X = mat2cell(X, 2*ones(Ltracers,1), nt);


%     timer = tic;
%     for i = 1:(nt-1)
%         if ~mod(i,100)
%            percentdone = 100 * i / (nt-1);
%            msg = sprintf('Percent done: %3.1f\n', percentdone);
%            fprintf(msg);
%         end
        
%         qq = fread(fidqh,2*Ntotal,'double');
%         qq = complex(qq(1:Ntotal), qq(Ntotal+1:2*Ntotal));
    
%         qh = reshape(qq,N,N);
%         psi_hat = InvBT.*qh;
    
%         v = real(ifft2(dX.*psi_hat));
%         u = real(ifft2(-dY.*psi_hat));
    
%         for n = 1:Ltracers
%             Xu = mod(X{n}(1,i), x(end));
%             Xv = mod(X{n}(2,i), x(end));
            
%             uxi = interp2(xx, yy, u, Xu, Xv);
%             vxi = interp2(xx, yy, v, Xu, Xv);
            
%             X{n}(1,i+1) = X{n}(1,i) + sigmax;
%             X{n}(2,i+1) = X{n}(2,i) + sigmax;
%         end
        
%     end
%     fprintf('\n tracer solve time %s\n', time2str(toc(timer)));
    
%     if fsavedata
%         casestr = sprintf('obss2_Ltracers%i_sigmax%.2g', Ltracers, sigmax);
%         save(fullfile(datafilename, strcat(casestr,'.mat')), 'X', 'sigmax', 'Ltracers', '-v7.3')
%     end
    
%     end