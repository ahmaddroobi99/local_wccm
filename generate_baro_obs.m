function X = generate_baro_obs(datafilename, Ltracers, sigmax, varargin)

p = inputParser;
addParameter(p, 'savedata', false, @islogical)

parse(p, varargin{:});
fsavedata = p.Results.savedata;


fnameqh = fullfile(datafilename, 'qh.bin');
fnamet = fullfile(datafilename, 't.bin');
dataext = load(fullfile(datafilename, 'diagnostics.mat'));

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


X = zeros(2*Ltracers, nt);
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

    parfor n = 1:Ltracers
        Xu = mod(X{n}(1,i), x(end));
        Xv = mod(X{n}(2,i), x(end));
        
        uxi = interp2(xx, yy, u, Xu, Xv);
        vxi = interp2(xx, yy, v, Xu, Xv);
        
        X{n}(1,i+1) = X{n}(1,i) + uxi*dt + sigmax*randn()*sqrt(dt);
        X{n}(2,i+1) = X{n}(2,i) + vxi*dt + sigmax*randn()*sqrt(dt);
    end
    
end
fprintf('\n tracer solve time %s\n', time2str(toc(timer)));

if fsavedata
    casestr = sprintf('obs_Ltracers%i_sigmax%.2g', Ltracers, sigmax);
    save(fullfile(datafilename, strcat(casestr,'.mat')), 'X', 'sigmax', 'Ltracers', '-v7.3')
end

end