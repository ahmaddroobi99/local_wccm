clear;

fsave = false;


casename = 'baroRK4_N128_beta0_d0.005_nu6.3e-25_forcing-constant_tau0-0_kf30_epsf0.05_sigf1_dt0.001_Nt200000_IC-supplied_set-10-29-23-085323';
datafolder = fullfile('C:\Users\ahmad\OneDrive\Desktop\Research\LTM\barotropic-tracers-master\barotropic-tracers-master\data',casename);
dataext = load(fullfile(datafolder, 'diagnostics.mat'));

Ltracers = 40;

obsdatafile = dir(fullfile(datafolder, sprintf('*Ltracers%i*.mat', Ltracers)));
if isempty(obsdatafile)
    error('obs data not found');
end

obsdata = load(fullfile(obsdatafile.folder, obsdatafile.name));
X = obsdata.X;

% X = generate_baro_obs(datafolder,Ltracers,sigmax,'savedata',true);

%% plot

figfolder = 'figs';
animdir = fullfile(figfolder,casename,'anim');

if fsave
    if ~isdir(animdir)
        mkdir(animdir)
    else
%         foldercontent = dir(animdir);
%         if length(foldercontent) > 2
%             prompt = 'Do you want to delete existing figure files? Y/N [Y]: ';
%             str = input(prompt,'s');
%             if isempty(str)
%                 str = 'Y';
%             end
%             if strcmp(str,'Y')
%                 rmdir(animdir); mkdir(animdir);
%             end
%         end
    end
end

fnameqh = fullfile(datafolder, 'qh.bin');
fnamet = fullfile(datafolder, 't.bin');

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


fplottracers = true;

figure;
XX = zeros(2,1);
ii = 0;
for i = 1:nt
    
    qq = fread(fidqh,2*Ntotal,'double');
    qq = complex(qq(1:Ntotal), qq(Ntotal+1:2*Ntotal));
    
    qh = reshape(qq,N,N);
    psi_hat = InvBT.*qh;
    psi = real(ifft2(psi_hat));
    
    if ~mod(i,10)
        ii = ii + 1;
        
        clf; hold on;
        
        pcolor(xx, yy, psi); view(2); view(2); shading interp;  box on; axis tight;
        set(gca,'XTick',[]);
        set(gca,'YTick',[]);
        colorbar('location','eastoutside', 'FontSize', 8);
        
%         caxis([-4 1.75])
        if fplottracers
            for n = 1:Ltracers
                XX(1) = mod(X{n}(1,i), x(end));
                XX(2) = mod(X{n}(2,i), x(end));
                plot(XX(1), XX(2), '.',  'MarkerSize', 14)
            end
        end
        xlim([0 x(end)]);
        ylim([0 x(end)]);
        set(gca,'Layer','top');
        
        title(sprintf('streamfunction - time = %.3g', t(i)));
        drawnow limitrate;
        if fsave
            export_fig(fullfile(animdir, sprintf('%07i',ii)), '-q101','-p.01','-m1','-png');
        end
    end
end


%%








