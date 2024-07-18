clear;

fsavefigs = true;

casename = 'baroRK4_N128_beta0.05_d0.005_nu2.5e-28_forcing-sin_tau0-0.04_dt0.02_Nt100_IC-supplied_set-11-13-23-134958';

datafolder = fullfile('data',casename);
dataext = load(fullfile('diag11.mat'));

figfolder = fullfile('figs_modes_full_data', casename);
if fsavefigs && ~isdir(figfolder); mkdir(figfolder); end


%%
fnameqh = fullfile('qh11.bin');
fnamet = fullfile('t11.bin');

fidqh = fopen(fnameqh,'r');
fidt = fopen(fnamet,'r');

t = fread(fidt,'double');
dt = t(2) - t(1);
nt = length(t);
N =  dataext.params.N;
Ntotal = N*N;

% 
% k = [0:N/2 -N/2+1:-1]';
% [KX, KY] = meshgrid(k,k);
% Knorm = hypot(KX,KY);
% Kinvnorm = 1./Knorm; Kinvnorm(1,1) = 0;
% dX = 1i*repmat(k', [N 1]);
% dY = 1i*repmat(k, [1 N]);
% Laplacian = dX.^2 + dY.^2;
% InvBT = 1./Laplacian; InvBT(1,1) = 0;
% 
% L = 2*pi;
% x = L/N*(0:N-1)';             % x coordinate
% y = L/N*(0:N-1)';             % y coordinate
% [xx,yy] = meshgrid(x,y);      % xy grid
% 
% 
% Ek = zeros(N/2+1,length(t));
% 
% qlmodes = zeros(4,nt); % save 1,0; 0,1; 1,1; -1,1 modes in order
% qlmodesa = zeros(4,nt); % save 1,0; 0,1; 1,1; -1,1 modes in order
% qlmodesb = zeros(4,nt);
% reversestr = '';
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % R&D 
% cases = loadcases();
% casename = cases(3);
% display(casename);
% datafolder = fullfile('data',casename);
% display(datafolder);
% 
% dataext = load(fullfile('diag11.mat'));
% display(dataext);
% 
% fnamet = fullfile('t11.bin'); fidt = fopen(fnamet,'r');
% display(fnamet);
% 
% t = fread(fidt,'double'); fclose(fidt);
% dt = t(2) - t(1);
% nt = length(t);
% Ntruth =  dataext.params.N;
% 
% 
% 
% 
% 
% 
% Ltracers = 20;
% sigmax = 0.20;
% obsdatafile = dir(fullfile(sprintf('obs_Ltracers%i_sigmax%.2g.mat', Ltracers, sigmax)));
% if isempty(obsdatafile)
%     error('obs data not found');
% end
% obsdata = load(fullfile(obsdatafile.folder,obsdatafile.name));
% 
% 
% pin.Tend = t(end);
% 
% Nf = 16;
% pin.Ltracers = Ltracers;
% pin.N = Nf; % must be even, this is the total grid points in each dimension, e.g  16    32    64  
% pin.X = obsdata.X;
% pin.E0 = 0.8;
% pin.damping_type = 'scale_dependent';
% pin.d = 0;
% pin.mu = 0.05;
% 
% pin.spectrum_type = 'third';
% pin.num_modes_flat = 0;
% 
% pin.sigmax = obsdata.sigmax;
% pin.Ltracers = obsdata.Ltracers;
% pin.dt = dt;
% pin.dtsave = dt*10;
% pin.fsavedata = true;
% pin.datafolder = casename;
% 
% 
% paramstrain = train_msm(datafolder,dt*20);
% params.train = 'msm';
% d = (paramstrain.d - 1i*paramstrain.omega);
% sigma = paramstrain.sigma;
% 
% k = [0:Ntruth/2 -Ntruth/2+1:-1]';
% [KX, KY] = meshgrid(k,k);
% 
% k = [0:pin.N/2 -pin.N/2+1:-1]';
% [KX, KY] = meshgrid(k,k);
% 
% pin.damping_type = 'scaledep';
% pin.spectrum_type = 'kolmo';
% 
% 
% [pin, m, Xt] = lagrangianfilter_baro(pin);
% qhc = coarsenqh(datafolder,pin.dtsave);
% 
% % display(m);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% for kk = 1:nt
%     if ~mod(kk,100)
%        percentdone = 100 * kk / nt;
%        msg = sprintf('Percent done: %3.1f', percentdone);
%        fprintf(strcat(reversestr, msg));
%        reversestr = repmat(sprintf('\b'), 1, length(msg));
%     end
% 
%     qq = fread(fidqh,2*Ntotal,'double');
%     % m=fread(m,2*Ntotal,'double');
%     qq = complex(qq(1:Ntotal), qq(Ntotal+1:2*Ntotal));
%     % mm = complex(m(1:Ntotal), m(Ntotal+1:2*Ntotal));
% 
%     qh = reshape(qq,N,N)/(N)^2;
%     % m = reshape(mm,N,N)/(N)^2;
% 
%     psi_hat = InvBT.*qh;
%     % psi_hat_2 = InvBT.*mm;
% 
%     vh = Knorm.*psi_hat;
%     % vh_2=Knorm.*psi_hat_2;
% 
%     qlmodes(1,kk) = vh(1,2); % 1,0
%     qlmodes(2,kk) = vh(2,1); % 0,1
%     qlmodes(3,kk) = vh(2,2); % 1,1
%     qlmodes(4,kk) = vh(2,N); % -1,1
% 
%     % conjugate modes (not needed)
%     qlmodesa(1,kk) = qh(1,N); % 1,0
%     qlmodesa(2,kk) = qh(N,1); % 0,1
%     qlmodesa(3,kk) = qh(N,N); % 1,1
%     qlmodesa(4,kk) = qh(N,2); % -1,1
% 
%     qlmodesb(1,kk) = m(1,N); % 1,0
%     qlmodesb(2,kk) = m(N,1); % 0,1
%     qlmodesb(3,kk) = m(N,N); % 1,1
%     qlmodesb(4,kk) = m(2,N); % -1,1
% 
% 
%     qlmodesc(1,kk) = qhc(1,N); % 1,0
%     qlmodesc(2,kk) = qhc(N,1); % 0,1
%     qlmodesc(3,kk) = qhc(N,N); % 1,1
%     qlmodesc(4,kk) = qhc(2,N); % -1,1
% 
% end
% %%
% handle = figure(1);
% clf;
% labels = ["(1,0)", "(0,1)", "(1,1)", "(-1,1)"];
% for i = 1:size(qlmodes,1)
%     subplot(4,1,i)
%     [Rmm,lags] = xcorr((qlmodes(i,1:20000).'- mean(qlmodes(i,1:20000),2).'),'coeff');
%     % [Rmm,lags] = xcorr((qlmodes(i,1:min(20000, end)).'- mean(qlmodes(i,1:min(20000, end)),2).'),'coeff');
% 
%     Rmm = Rmm(lags>0);
%     lags = lags(lags>0);
%     plot(lags*dt,real(Rmm)); hold on;
%     % plot(lags*dt,imag(Rmm))
%     % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%     % 
%     %    [Rmm,lags] = xcorr((qlmodesa(i,1:20000).'- mean(qlmodesa(i,1:20000),2).'),'coeff');
%     % % [Rmm,lags] = xcorr((qlmodes(i,1:min(20000, end)).'- mean(qlmodes(i,1:min(20000, end)),2).'),'coeff');
%     % 
%     % Rmm = Rmm(lags>0);
%     % lags = lags(lags>0);
%     % plot(lags*dt,real(Rmm)); hold on;
%     % % plot(lags*dt,imag(Rmm))
% 
% 
% 
%     if i == 1 || i == 2
%     xlim([0 5e3*dt])
%     else
%     xlim([0 1e4*dt]); 
%     end
%     title(sprintf('mode %s',labels(i)));
%     xlabel('t');
% end
% set(gcf,'Position',[867 204 348 709]);
% % % % % % % % % % % % % % % % 
% if fsavefigs
%     resizefonts(handle);
%     export_fig(fullfile(figfolder, sprintf('autocorr2')), '-q101','-p.01','-m1','-pdf');
%     savefig(handle, fullfile(figfolder, sprintf('autocorr2.fig')))
% end
% handle2 = figure(2);
% 
% labels = ["(1,0)", "(0,1)", "(1,1)", "(-1,1)"];
% for i = 1:size(qlmodesa,1)
%     subplot(4,1,i)
%     % [Rmm,lags] = xcorr((qlmodes(i,1:20000).'- mean(qlmodes(i,1:20000),2).'),'coeff');
%     % % [Rmm,lags] = xcorr((qlmodes(i,1:min(20000, end)).'- mean(qlmodes(i,1:min(20000, end)),2).'),'coeff');
%     % 
%     % Rmm = Rmm(lags>0);
%     % lags = lags(lags>0);
%     % plot(lags*dt,real(Rmm)); hold on;
%     % % plot(lags*dt,imag(Rmm))
%     % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% 
%        [Rmm,lags] = xcorr((qlmodesa(i,1:20000).'- mean(qlmodesa(i,1:20000),2).'),'coeff');
%     % [Rmm,lags] = xcorr((qlmodes(i,1:min(20000, end)).'- mean(qlmodes(i,1:min(20000, end)),2).'),'coeff');
% 
%     Rmm = Rmm(lags>0);
%     lags = lags(lags>0);
%     plot(lags*dt,real(Rmm)); hold on;
%     % plot(lags*dt,imag(Rmm))
% 
% 
% 
%     if i == 1 || i == 2
%     xlim([0 5e3*dt])
%     else
%     xlim([0 1e4*dt]); 
%     end
%     title(sprintf('mode %s',labels(i)));
%     xlabel('t');
% end
% % set(gcf,'Position',[867 204 348 709]);
% 
% if fsavefigs
%     resizefonts(handle2);
%     export_fig(fullfile(figfolder, sprintf('autocorr')), '-q101','-p.01','-m1','-pdf');
%     savefig(handle2, fullfile(figfolder, sprintf('autocorr.fig')))
% end
% 
% 
% 
% 
% handle3 = figure(3);
% clf;
% labels = ["(1,0)", "(0,1)", "(1,1)", "(-1,1)"];
% for i = 1:size(qlmodesb,1)
%     subplot(4,1,i)
%     [Rmm,lags] = xcorr((qlmodesb(i,1:20000).'- mean(qlmodesb(i,1:20000),2).'),'coeff');
%     % [Rmm,lags] = xcorr((qlmodes(i,1:min(20000, end)).'- mean(qlmodes(i,1:min(20000, end)),2).'),'coeff');
% 
%     Rmm = Rmm(lags>0);
%     lags = lags(lags>0);
%     plot(lags*dt,real(Rmm)); hold on;
%     % plot(lags*dt,imag(Rmm))
%     % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%     % 
%     %    [Rmm,lags] = xcorr((qlmodesa(i,1:20000).'- mean(qlmodesa(i,1:20000),2).'),'coeff');
%     % % [Rmm,lags] = xcorr((qlmodes(i,1:min(20000, end)).'- mean(qlmodes(i,1:min(20000, end)),2).'),'coeff');
%     % 
%     % Rmm = Rmm(lags>0);
%     % lags = lags(lags>0);
%     % plot(lags*dt,real(Rmm)); hold on;
%     % % plot(lags*dt,imag(Rmm))
% 
% 
% 
%     if i == 1 || i == 2
%     xlim([0 5e3*dt])
%     else
%     xlim([0 1e4*dt]); 
%     end
%     title(sprintf('mode %s',labels(i)));
%     xlabel('t');
% end
% % set(gcf,'Position',[867 204 348 709]);
% 
% if fsavefigs
%     resizefonts(handle3);
%     export_fig(fullfile(figfolder, sprintf('autocorr3')), '-q101','-p.01','-m1','-pdf');
%     savefig(handle3, fullfile(figfolder, sprintf('autocorr3.fig')))
% end
% 
% 
% 
% 
% 
% handle4 = figure(4);
% clf;
% labels = ["(1,0)", "(0,1)", "(1,1)", "(-1,1)"];
% for i = 1:size(qlmodesb,1)
%     subplot(4,1,i)
%     [Rmm,lags] = xcorr((qlmodesc(i,1:20000).'- mean(qlmodesc(i,1:20000),2).'),'coeff');
%     % [Rmm,lags] = xcorr((qlmodes(i,1:min(20000, end)).'- mean(qlmodes(i,1:min(20000, end)),2).'),'coeff');
% 
%     Rmm = Rmm(lags>0);
%     lags = lags(lags>0);
%     plot(lags*dt,real(Rmm)); hold on;
%     % plot(lags*dt,imag(Rmm))
%     % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%     % 
%     %    [Rmm,lags] = xcorr((qlmodesa(i,1:20000).'- mean(qlmodesa(i,1:20000),2).'),'coeff');
%     % % [Rmm,lags] = xcorr((qlmodes(i,1:min(20000, end)).'- mean(qlmodes(i,1:min(20000, end)),2).'),'coeff');
%     % 
%     % Rmm = Rmm(lags>0);
%     % lags = lags(lags>0);
%     % plot(lags*dt,real(Rmm)); hold on;
%     % % plot(lags*dt,imag(Rmm))
% 
% 
% 
%     if i == 1 || i == 2
%     xlim([0 5e3*dt])
%     else
%     xlim([0 1e4*dt]); 
%     end
%     title(sprintf('mode %s',labels(i)));
%     xlabel('t');
% end
% % set(gcf,'Position',[867 204 348 709]);
% 
% if fsavefigs
%     resizefonts(handle4);
%     export_fig(fullfile(figfolder, sprintf('autocorr4')), '-q101','-p.01','-m1','-pdf');
%     savefig(handle4, fullfile(figfolder, sprintf('autocorr4.fig')))
% end





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


Ek = zeros(N/2+1,length(t));

qlmodes = zeros(4,nt); % save 1,0; 0,1; 1,1; -1,1 modes in order
% qlmodesa = zeros(4,nt); % save 1,0; 0,1; 1,1; -1,1 modes in order
reversestr = '';
for kk = 1:nt
    if ~mod(kk,100)
       percentdone = 100 * kk / nt;
       msg = sprintf('Percent done: %3.1f', percentdone);
       fprintf(strcat(reversestr, msg));
       reversestr = repmat(sprintf('\b'), 1, length(msg));
    end
    
    qq = fread(fidqh,2*Ntotal,'double');
    qq = complex(qq(1:Ntotal), qq(Ntotal+1:2*Ntotal));

    qh = reshape(qq,N,N)/(N)^2;
    psi_hat = InvBT.*qh;
    vh = Knorm.*psi_hat;
    
    qlmodes(1,kk) = vh(1,2); % 1,0
    qlmodes(2,kk) = vh(2,1); % 0,1
    qlmodes(3,kk) = vh(2,2); % 1,1
    qlmodes(4,kk) = vh(2,N); % -1,1
    
    % conjugate modes (not needed)
%     qlmodesa(1,kk) = qh(1,N); % 1,0
%     qlmodesa(2,kk) = qh(N,1); % 0,1
%     qlmodesa(3,kk) = qh(N,N); % 1,1
%     qlmodesa(4,kk) = qh(N,2); % -1,1

end
%%
handle = figure(1);
clf;
labels = ["(1,0)", "(0,1)", "(1,1)", "(-1,1)"];
for i = 1:size(qlmodes,1)
    subplot(4,1,i)
    [Rmm,lags] = xcorr((qlmodes(i,1:20000).'- mean(qlmodes(i,1:20000),2).'),'coeff');
    Rmm = Rmm(lags>0);
    lags = lags(lags>0);
    plot(lags*dt,real(Rmm)); hold on;
%     plot(lags*dt,imag(Rmm))
    if i == 1 || i == 2
    xlim([0 5e4*dt])
    else
    xlim([0 1e4*dt]); 
    end
    title(sprintf('mode %s',labels(i)));
    xlabel('t');
end
set(gcf,'Position',[867 204 348 709]);


if fsavefigs
    resizefonts(handle);
    export_fig(fullfile(figfolder, sprintf('autocorr')), '-q101','-p.01','-m1','-pdf');
    savefig(handle, fullfile(figfolder, sprintf('autocorr.fig')))
end