% 
% cases = 's';
% 
% fsavefigs = true;
% 
% for casename = cases
% 
% 
% datafolder = fullfile('data',char(casename));
% figfolder = fullfile('figs_spectra_xx', char(casename));
% if fsavefigs && ~isdir(figfolder); mkdir(figfolder); end
% 
% dataext = load(fullfile('diag11.mat'));
% 
% qhc = coarsenqh(datafolder, 0.02); % load data, second argument is junk
% 
% Ntruth = sqrt(size(qhc,1));
% k = [0:Ntruth/2 -Ntruth/2+1:-1]';
% [KX, KY] = meshgrid(k,k);
% Knorm = hypot(KX,KY);
% Knorminv = 1./Knorm; Knorminv(1,1) = 0;
% 
% kp = Knorm(1,1:Ntruth/2+1);
% 
% 
% psi_hat_qhc = -Knorminv(:).*qhc/Ntruth^2;
% qhcm = abs(mean(psi_hat_qhc,2)).^2;
% qhcv = var(psi_hat_qhc,0,2);
% 
% qhcm = reshape(qhcm,Ntruth,Ntruth);
% qhcv = reshape(qhcv,Ntruth,Ntruth);
% 
% ev = zeros(Ntruth/2+1,1);
% em = zeros(Ntruth/2+1,1);
% for jj = 1:Ntruth
%     for ii = 1:Ntruth
%         k = sqrt(KX(ii,jj).^2 + KY(ii,jj).^2);
%         if ceil(k) <= Ntruth/2
%             r = k-floor(k);
%             % plus 1 because of matlab indexing
%             ev(floor(k)+1) = ev(floor(k)+1) + (1-r)*abs(qhcv(ii,jj));
%             ev(ceil(k)+1) = ev(ceil(k)+1) + r*abs(qhcv(ii,jj));
% 
%             em(floor(k)+1) = em(floor(k)+1) + (1-r)*abs(qhcm(ii,jj));
%             em(ceil(k)+1) = em(ceil(k)+1) + r*abs(qhcm(ii,jj));
%         end
%     end
% end
% ev = 0.5*ev;
% em = 0.5*em;
% 
% 
% %%
% handle = figure;
% loglog(kp,em+ev,'LineWidth',1.25);
% hold on;
% 
% loglog(kp,ev); 
% loglog(kp,em);
% legend('total','variance','mean');
% 
% 
% if strcmp(dataext.params.forcingtype,'constant')
%     kkp = (1:35)';
% 
%     alpha = -5/3;
% else
%         kkp = (1:45)';
% 
%     alpha = -3;
% end
% [~,~,ikp] = intersect(kkp,kp);
% amin = fminsearch(@(a) mean(log(abs(a*kkp.^alpha - ev(ikp)))), 1);
% 
% loglog(kkp,amin*kkp.^alpha,'k--');
% text(2,amin*kkp(1).^alpha,strcat('k^{',strtrim(rats(alpha)),'}'))
% 
% 
% 
% xlabel('wavenumber');
% ylabel('energy')
% 
% drawnow; 
% if fsavefigs
%     resizefonts(handle);
%     export_fig(fullfile(figfolder, sprintf('spectrum')), '-q101','-p.01','-m1','-pdf');
%     savefig(handle, fullfile(figfolder, sprintf('spectrum.fig')))
% end
% 
% 
% 
% end

cases = 's'; % Assuming 'cases' is a character array, not a cell array
fsavefigs = true;

for casename = cases
    
    datafolder = fullfile('data', casename);
    figfolder = fullfile('figs_spectra_xx', casename);
    if fsavefigs && ~isdir(figfolder)
        mkdir(figfolder);
    end
    
    dataext = load(fullfile('diag11.mat'));
    
    qhc = coarsenqh(datafolder, 0.02); % load data, second argument is junk
    
    Ntruth = sqrt(size(qhc, 1));
    k = [0:Ntruth/2, -Ntruth/2+1:-1]';
    [KX, KY] = meshgrid(k, k);
    Knorm = hypot(KX, KY);
    Knorminv = 1 ./ Knorm;
    Knorminv(1, 1) = 0;
    
    kp = Knorm(1, 1:Ntruth/2+1);
    
    psi_hat_qhc = -Knorminv(:) .* qhc / Ntruth^2;
    qhcm = abs(mean(psi_hat_qhc, 2)).^2;
    qhcv = var(psi_hat_qhc, 0, 2);
    
    qhcm = reshape(qhcm, Ntruth, Ntruth);
    qhcv = reshape(qhcv, Ntruth, Ntruth);
    
    ev = zeros(Ntruth/2+1, 1);
    em = zeros(Ntruth/2+1, 1);
    for jj = 1:Ntruth
        for ii = 1:Ntruth
            k = sqrt(KX(ii, jj).^2 + KY(ii, jj).^2);
            if ceil(k) <= Ntruth/2
                r = k-floor(k);
                % plus 1 because of MATLAB indexing
                ev(floor(k)+1) = ev(floor(k)+1) + (1-r) * abs(qhcv(ii, jj));
                ev(ceil(k)+1) = ev(ceil(k)+1) + r * abs(qhcv(ii, jj));
                
                em(floor(k)+1) = em(floor(k)+1) + (1-r) * abs(qhcm(ii, jj));
                em(ceil(k)+1) = em(ceil(k)+1) + r * abs(qhcm(ii, jj));
            end
        end
    end
    ev = 0.5 * ev;
    em = 0.5 * em;
    
    % Divide data into three segments
    seg_length = floor(length(kp) / 3);
    kp_segments = {kp(1:seg_length), kp(1:2*seg_length), kp(1:end)};
    em_segments = {em(1:seg_length), em(1:2*seg_length), em(1:end)};
    ev_segments = {ev(1:seg_length), ev(1:2*seg_length), ev(1:end)};
    
    for seg_idx = 1:3
        handle = figure;
        loglog(kp_segments{seg_idx}, em_segments{seg_idx} + ev_segments{seg_idx}, 'LineWidth', 1.25);
        hold on;
        loglog(kp_segments{seg_idx}, ev_segments{seg_idx}); 
        loglog(kp_segments{seg_idx}, em_segments{seg_idx});
        legend('total', 'variance', 'mean');
        
        if strcmp(dataext.params.forcingtype, 'constant')
            kkp = (1:35)';
            alpha = -5/3;
        else
            kkp = (1:45)';
            alpha = -3;
        end
        [~, ~, ikp] = intersect(kkp, 30);
        kkp_interp = kkp(ikp); % Interpolate kkp to match dimensions with ikp
        ev_interp = interp1(kp_segments{seg_idx}, ev_segments{seg_idx}, kkp_interp); % Interpolate ev_segments to match dimensions with kkp_interp
        amin = fminsearch(@(a) mean(log(abs(a * kkp_interp.^alpha - ev_interp))), 1);
                
        loglog(kkp, amin * kkp.^alpha, 'k--');
        text(2, amin * kkp(1).^alpha, strcat('k^{', strtrim(rats(alpha)), '}'));
        
        xlabel('wavenumber');
        ylabel('energy');
        
        drawnow; 
        if fsavefigs
            resizefonts(handle);
            export_fig(fullfile(figfolder, sprintf('spectrum_segment%d', seg_idx)), '-q101', '-p.01', '-m1', '-pdf');
            savefig(handle, fullfile(figfolder, sprintf('spectrum_segment%d.fig', seg_idx)))
        end
        
        % Pause to allow observation of the figure before proceeding
        pause;
    end
end

