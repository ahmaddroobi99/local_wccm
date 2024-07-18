display("Let's code the Time Series here!");

clear;
% % % % some changes to test git in arc
fsavefigs = true;

casename = 'baroRK4_N128_beta0.05_d0.005_nu2.5e-28_forcing-sin_tau0-0.04_dt0.02_Nt100_IC-supplied_set-11-13-23-134958';

datafolder = fullfile('data', casename);
dataext = load(fullfile('diag11.mat'));

figfolder = fullfile('figs_time_series_V3', casename);
if fsavefigs && ~isdir(figfolder); mkdir(figfolder); end

%%
fnameqh = fullfile('qh11.bin');
fnamet = fullfile('t11.bin');

fidqh = fopen(fnameqh, 'r');
display(size(fidqh));
fidt = fopen(fnamet, 'r');

t = fread(fidt, 'double');
dt = t(2) - t(1);
nt = length(t);
display(size(t));

N = dataext.params.N;
Ntruth = dataext.params.N;

Ntotal = N * N;

k = [0:N/2 -N/2+1:-1]';
display(size(k));

[KX, KY] = meshgrid(k, k);

Knorm = hypot(KX, KY);
display(size(Knorm));

Kinvnorm = 1./Knorm;
display(size(Kinvnorm));

Kinvnorm(1,1) = 0;
display(size(Kinvnorm));

dX = 1i * repmat(k', [N 1]);
display(size(dX));

dY = 1i * repmat(k, [1 N]);
Laplacian = dX.^2 + dY.^2;
display(size(Laplacian));

InvBT = 1./Laplacian;
display(size(InvBT));

InvBT(1,1) = 0;

L = 2 * pi;
x = L/N * (0:N-1)';             % x coordinate
y = L/N * (0:N-1)';             % y coordinate
[xx, yy] = meshgrid(x, y);      % xy grid

Ek = zeros(N/2+1, length(t));

qlmodes = zeros(4, nt); % save 1,0; 0,1; 1,1; -1,1 modes in order
mm = zeros(4, nt); % save 1,0; 0,1; 1,1; -1,1 modes in order
reversestr = '';
for kk = 1:nt
    % if ~mod(kk, 100)
    %     percentdone = 100 * kk / nt;
    %     msg = sprintf('Percent done: %3.1f', percentdone);
    %     fprintf(strcat(reversestr, msg));
    %     reversestr = repmat(sprintf('\b'), 1, length(msg));
    % end
    
    qq = fread(fidqh, 2 * Ntotal, 'double');
    display(size(qq));

    qq = complex(qq(1:Ntotal), qq(Ntotal+1:2*Ntotal));
    display(size(qq));

    qh = reshape(qq, N, N) / (N)^2;
        display(size(qh));

    psi_hat = InvBT .* qh;
     display(size(psi_hat));

    vh = Knorm .* psi_hat;
         display(size(vh));

    qlmodes(1,kk) = vh(1,2); % 1,0
    qlmodes(2,kk) = vh(2,1); % 0,1
             display(size(    qlmodes(2,kk)));

    qlmodes(3,kk) = vh(2,2); % 1,1
    qlmodes(4,kk) = vh(2,N); % -1,1
         display(size(    qlmodes));

    % conjugate modes (not needed)
%     qlmodesa(1,kk) = qh(1,N); % 1,0
%     qlmodesa(2,kk) = qh(N,1); % 0,1
%     qlmodesa(3,kk) = qh(N,N); % 1,1
%     qlmodesa(4,kk) = qh(N,2); % -1,1

end
%%
% handle = figure(1);
% clf;
% labels = ["(1,0)", "(0,1)", "(1,1)", "(-1,1)"];
% for i = 1:size(qlmodes,1)
%     subplot(4,1,i)
%     % [Rmm,lags] = xcorr((qlmodes(i,1:20000).'- mean(qlmodes(i,1:20000),2).'),'coeff');
%     [Rmm,lags] = xcorr((qlmodes(i,1:min(20000, end)).'- mean(qlmodes(i,1:min(20000, end)),2).'),'coeff');
% 
%     Rmm = Rmm(lags>0);
%     lags = lags(lags>0);
%     plot(lags*dt,real(Rmm)); hold on;
% %     plot(lags*dt,imag(Rmm))
%     if i == 1 || i == 2
%     xlim([0 5e3*dt])
%     else
%     xlim([0 1e4*dt]); 
%     end
%     title(sprintf('mode %s',labels(i)));
%     xlabel('t');
% end
% set(gcf,'Position',[867 204 348 709]);

% 
% if fsavefigs
%     resizefonts(handle);
%     export_fig(fullfile(figfolder, sprintf('autocorr')), '-q101','-p.01','-m1','-pdf');
%     savefig(handle, fullfile(figfolder, sprintf('autocorr.fig')))
% end


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% R&D

%%%%%%%%%%%%%%%%%
% observations
%%%%%%%%%%%%%%%%%
Ltracers = 20;
sigmax = 0.20;
obsdatafile = dir(fullfile(sprintf('obs_Ltracers%i_sigmax%.2g.mat', Ltracers, sigmax)));
if isempty(obsdatafile)
    error('obs data not found');
end
obsdata = load(fullfile(obsdatafile.folder,obsdatafile.name));

pin.Tend = t(end);

Nf = 16;
pin.Ltracers = Ltracers;
pin.N = Nf; % must be even, this is the total grid points in each dimension, e.g  16    32    64  
pin.X = obsdata.X;
pin.E0 = 0.8;
pin.damping_type = 'scaledep';
pin.d = 0;
pin.mu = 0.05;

pin.spectrum_type = 'third';
pin.num_modes_flat = 0;

pin.sigmax = obsdata.sigmax;
pin.Ltracers = obsdata.Ltracers;
pin.dt = dt;
pin.dtsave = dt * 10;
pin.fsavedata = true;
pin.datafolder = casename;


paramstrain = train_msm(datafolder,dt*20);
params.train = 'msm';
d = (paramstrain.d - 1i * paramstrain.omega);
sigma = paramstrain.sigma;

k = [0:Ntruth/2 -Ntruth/2+1:-1]';
[KX, KY] = meshgrid(k, k);

k = [0:pin.N/2 -pin.N/2+1:-1]';
[KX, KY] = meshgrid(k, k);

% pin.damping_type = 'scaledep';
% pin.spectrum_type = 'kolmo';


[pin, m, Xt] = lagrangianfilter_baro(pin);


display(size(m));
% % % % % % % % % % % % % % % % % % % % % % % 
for kk = 1:nt
    if ~mod(kk, 100)
       percentdone = 100 * kk / nt;
       % msg = sprintf('Percent done: %3.1f', percentdone);
       % fprintf(strcat(reversestr, msg));
       % reversestr = repmat(sprintf('\b'), 1, length(msg));
    end
    
    % Reshape m to NxN
% mmm = reshape(m(:, kk), [], N);
    % psi_hat = InvBT.*m;
    % m = Knorm.*psi_hat;
    mmm=m;
    % display(mmm);
    % 
    % mm(1, kk) = mmm(1, 2); % 1,0
    % mm(2, kk) = mmm(2, 1); % 0,1
    % mm(3, kk) = mmm(2, 2); % 1,1
    % mm(4, kk) = mmm(2, N); % -1,1

end

% % % % % % % % % % % % % % % % % % % % % % % 

handle = figure(1);
clf;
labels = ["(1,0)", "(0,1)", "(1,1)", "(-1,1)"];

time = dt * (0:size(qlmodes, 2) - 1); % Time vector

for i = 1:size(mmm, 1)
    subplot(4, 1, i)
    plot(time, real(qlmodes(i, :)), 'DisplayName', 'True'); % Plot real part of the mode (True)
    hold on;
    plot(time, (mmm(i, :)), 'DisplayName', 'Recovery'); % Plot real part of the mode (Recovery)
    title(sprintf('mode %s', labels(i)));
    xlabel('Time');
    ylabel('Amplitude');
    legend('show');
end

set(gcf, 'Position', [867 204 348 709]);

if fsavefigs
    resizefonts(handle);
    export_fig(fullfile(figfolder, sprintf('time_series')), '-q101', '-p.01', '-m1', '-pdf');
    savefig(handle, fullfile(figfolder, sprintf('time_series.fig')))
end
