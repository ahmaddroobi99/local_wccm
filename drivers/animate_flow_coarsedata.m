clear;

fsave = false;

cases = loadcases();

casename = cases(2);
datafolder = fullfile('C:\Users\ahmad\OneDrive\Desktop\Research\LTM\barotropic-tracers-master\barotropic-tracers-master\data',casename);
dataext = load(fullfile(datafolder, 'diagnostics.mat'));

searchfilename = dir(fullfile(datafolder, 'coarseqh_dt*.mat')); % arbitrarily load regardless of dtsave size
if isempty(searchfilename)
    error('not found');
end
    
data = load(fullfile(datafolder, searchfilename.name));


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

t =  data.tsave;
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


figure;
ii = 0;
for i = 1:nt
     
    qh = reshape(data.qhc(:,i),N,N);
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








