% 
% cases = 'spectra';
% 
% fsavefigs = true;
% 
% for casename = cases
% 
% 
% datafolder = fullfile('data',char(casename));
% figfolder = fullfile('figs_spectra', char(casename));
% if fsavefigs && ~isdir(figfolder); mkdir(figfolder); end
% 
% dataext = load(fullfile('diag11.mat'));
% 
% % qhc = coarsenqh(datafolder, 0.02); % load data, second argument is junk
% % 
% % Ntruth = sqrt(size(qhc,1));
% % k = [0:Ntruth/2 -Ntruth/2+1:-1]';
% % [KX, KY] = meshgrid(k,k);
% % Knorm = hypot(KX,KY);
% % Knorminv = 1./Knorm; Knorminv(1,1) = 0;
% % 
% % kp = Knorm(1,1:Ntruth/2+1);
% % 
% % 
% % psi_hat_qhc = -Knorminv(:).*qhc/Ntruth^2;
% % qhcm = abs(mean(psi_hat_qhc,2)).^2;
% % qhcv = var(psi_hat_qhc,0,2);
% % 
% % qhcm = reshape(qhcm,Ntruth,Ntruth);
% % qhcv = reshape(qhcv,Ntruth,Ntruth);
% % 
% % ev = zeros(Ntruth/2+1,1);
% % em = zeros(Ntruth/2+1,1);
% % for jj = 1:Ntruth
% %     for ii = 1:Ntruth
% %         k = sqrt(KX(ii,jj).^2 + KY(ii,jj).^2);
% %         if ceil(k) <= Ntruth/2
% %             r = k-floor(k);
% %             % plus 1 because of matlab indexing
% %             ev(floor(k)+1) = ev(floor(k)+1) + (1-r)*abs(qhcv(ii,jj));
% %             ev(ceil(k)+1) = ev(ceil(k)+1) + r*abs(qhcv(ii,jj));
% % 
% %             em(floor(k)+1) = em(floor(k)+1) + (1-r)*abs(qhcm(ii,jj));
% %             em(ceil(k)+1) = em(ceil(k)+1) + r*abs(qhcm(ii,jj));
% %         end
% %     end
% % end
% % ev = 0.5*ev;
% % em = 0.5*em;
% % 
% % 
% % %%
% % handle = figure;
% % loglog(kp,em+ev,'LineWidth',1.25);
% % hold on;
% % 
% % loglog(kp,ev); 
% % loglog(kp,em);
% % legend('total','variance','mean');
% % 
% % 
% % if strcmp(dataext.params.forcingtype,'constant')
% %     kkp = (1:35)';
% % 
% %     alpha = -5/3;
% % else
% %         kkp = (1:45)';
% % 
% %     alpha = -3;
% % end
% % [~,~,ikp] = intersect(kkp,kp);
% % amin = fminsearch(@(a) mean(log(abs(a*kkp.^alpha - ev(ikp)))), 1);
% % 
% % loglog(kkp,amin*kkp.^alpha,'k--');
% % text(2,amin*kkp(1).^alpha,strcat('k^{',strtrim(rats(alpha)),'}'))
% % 
% % 
% % 
% % xlabel('wavenumber');
% % ylabel('energy')
% % 
% % drawnow; 
% % if fsavefigs
% %     resizefonts(handle);
% %     export_fig(fullfile(figfolder, sprintf('spectrum')), '-q101','-p.01','-m1','-pdf');
% %     savefig(handle, fullfile(figfolder, sprintf('spectrum.fig')))
% % end
% % 
% % 
% % 
% % end
% 
% 
% 
% 
% qhc = coarsenqh(datafolder, 0); % load data, second argument is junk
% 
%     Ntruth = sqrt(size(qhc, 1));
%     k = [0:Ntruth/2 -Ntruth/2+1:-1]';
%     [KX, KY] = meshgrid(k, k);
%     Knorm = hypot(KX, KY);
%     Knorminv = 1 ./ Knorm;
%     Knorminv(1, 1) = 0;
% 
%     kp = Knorm(1, 1:Ntruth/2+1);
% 
%     num_segments = 10; % Define the number of segments
%     step_size = floor(size(qhc, 1) / num_segments); % Calculate the step size for segments
% 
%     em = zeros(Ntruth/2+1, num_segments);
%     ev = zeros(Ntruth/2+1, num_segments);
% 
%     for seg = 1:num_segments
%         data_segment = qhc(1:step_size*seg, :);
% 
%         psi_hat_qhc_segment = -Knorminv(1:step_size*seg) .* data_segment / (Ntruth*seg)^2;
%         qhcm_segment = abs(mean(psi_hat_qhc_segment, 2)).^2;
%         qhcv_segment = var(psi_hat_qhc_segment, 0, 2);
% 
%         % Reshape qhcm_segment and qhcv_segment based on the segment size
%         segment_size = sqrt(numel(qhcm_segment));
%         segment_size = ceil(segment_size);
% 
%         num_elements = segment_size * segment_size;
%         if num_elements > numel(qhcm_segment)
%             qhcm_segment = [qhcm_segment; zeros(num_elements - numel(qhcm_segment), 1)]; % Pad with zeros
%             qhcv_segment = [qhcv_segment; zeros(num_elements - numel(qhcv_segment), 1)]; % Pad with zeros
%         end
% 
%         qhcm_segment = reshape(qhcm_segment(1:num_elements), segment_size, segment_size);
%         qhcv_segment = reshape(qhcv_segment(1:num_elements), segment_size, segment_size);
% 
%         for jj = 1:segment_size
%             for ii = 1:segment_size
%                 k = sqrt(KX(ii, jj).^2 + KY(ii, jj).^2);
%                 if ceil(k) <= Ntruth/2
%                     r = k - floor(k);
% 
%                     idx_floor = min(max(1, floor(k)) + 1, size(ev, 1));
%                     idx_ceil = min(max(1, ceil(k)) + 1, size(ev, 1));
% 
%                     if idx_floor >= 1 && idx_floor <= size(ev, 1) && idx_ceil >= 1 && idx_ceil <= size(ev, 1)
%                         ev(idx_floor, seg) = ev(idx_floor, seg) + (1 - r) * abs(qhcv_segment(ii, jj));
%                         ev(idx_ceil, seg) = ev(idx_ceil, seg) + r * abs(qhcv_segment(ii, jj));
% 
%                         em(idx_floor, seg) = em(idx_floor, seg) + (1 - r) * abs(qhcm_segment(ii, jj));
%                         em(idx_ceil, seg) = em(idx_ceil, seg) + r * abs(qhcm_segment(ii, jj));
%                     end
%                 end
%             end
%         end
%         ev(:, seg) = 0.5 * ev(:, seg);
%         em(:, seg) = 0.5 * em(:, seg);
%     end
% 
%     % Calculate average of segments
%     ev_avg = mean(ev, 2);
%     em_avg = mean(em, 2);
% 
%     handle = figure;
%     loglog(kp, em_avg + ev_avg, 'LineWidth', 1.25);
%     hold on;
% 
%     loglog(kp, ev_avg);
%     loglog(kp, em_avg);
%     legend('total', 'variance', 'mean');
% 
%   % % % % % % % % % % % % % % % % 
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
% end


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% 
% cases = loadcases();
% 
% fsavefigs = true;
% 
% for casename = cases
% 
% 
% datafolder = fullfile('data',char(casename));
% figfolder = fullfile('figs', char(casename));
% if fsavefigs && ~isdir(figfolder); mkdir(figfolder); end
% 
% dataext = load(fullfile(datafolder, 'diagnostics.mat'));
% 
% qhc = coarsenqh(datafolder, 0); % load data, second argument is junk
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
% 
% end

% cases = loadcases(3);
% fsavefigs = true;
% 
% % for casename = cases
% 
%     datafolder = fullfile('data', '');
%     figfolder = fullfile('figs_spec__Ssss', '');
%     if fsavefigs && ~isdir(figfolder)
%         mkdir(figfolder);
%     % end
% 
%     dataext = load(fullfile('diag11.mat'));
% 
%     qhc = coarsenqh(datafolder, 0); % load data, second argument is junk
% 
%     Ntruth = sqrt(size(qhc, 1));
%     k = [0:Ntruth/2 -Ntruth/2+1:-1]';
%     [KX, KY] = meshgrid(k, k);
%     Knorm = hypot(KX, KY);
%     Knorminv = 1 ./ Knorm;
%     Knorminv(1, 1) = 0;
% 
%     kp = Knorm(1, 1:Ntruth/2+1);
% 
%     num_segments = 10;
%     step_size = floor(size(qhc, 1) / num_segments);
% 
%     em = zeros(Ntruth/2+1, num_segments);
%     ev = zeros(Ntruth/2+1, num_segments);
% 
%     for seg = 1:num_segments
%         data_segment = qhc(1:step_size*seg, :);
% 
%         psi_hat_qhc_segment = -Knorminv(1:step_size*seg) .* data_segment / (Ntruth*seg)^2;
%         qhcm_segment = abs(mean(psi_hat_qhc_segment, 2)).^2;
%         qhcv_segment = var(psi_hat_qhc_segment, 0, 2);
% 
%         % Reshape qhcm_segment and qhcv_segment based on the segment size
%         segment_size = sqrt(numel(qhcm_segment));
%         segment_size = ceil(segment_size);
% 
%         num_elements = segment_size * segment_size;
%         if num_elements > numel(qhcm_segment)
%             qhcm_segment = [qhcm_segment; zeros(num_elements - numel(qhcm_segment), 1)]; % Pad with zeros
%             qhcv_segment = [qhcv_segment; zeros(num_elements - numel(qhcv_segment), 1)]; % Pad with zeros
%         end
% 
%         qhcm_segment = reshape(qhcm_segment(1:num_elements), segment_size, segment_size);
%         qhcv_segment = reshape(qhcv_segment(1:num_elements), segment_size, segment_size);
% 
%         for jj = 1:segment_size
%             for ii = 1:segment_size
%                 k = sqrt(KX(ii, jj).^2 + KY(ii, jj).^2);
%                 if ceil(k) <= Ntruth/2
%                     r = k - floor(k);
% 
%                     idx_floor = min(max(1, floor(k)) + 1, size(ev, 1));
%                     idx_ceil = min(max(1, ceil(k)) + 1, size(ev, 1));
% 
%                     if idx_floor >= 1 && idx_floor <= size(ev, 1) && idx_ceil >= 1 && idx_ceil <= size(ev, 1)
%                         ev(idx_floor, seg) = ev(idx_floor, seg) + (1 - r) * abs(qhcv_segment(ii, jj));
%                         ev(idx_ceil, seg) = ev(idx_ceil, seg) + r * abs(qhcv_segment(ii, jj));
% 
%                         em(idx_floor, seg) = em(idx_floor, seg) + (1 - r) * abs(qhcm_segment(ii, jj));
%                         em(idx_ceil, seg) = em(idx_ceil, seg) + r * abs(qhcm_segment(ii, jj));
%                     end
%                 end
%             end
%         end
%         ev(:, seg) = 0.5 * ev(:, seg);
%         em(:, seg) = 0.5 * em(:, seg);
%     end
% 
%     % Calculate average of segments
%     ev_avg = mean(ev, 2);
%     em_avg = mean(em, 2);
% 
%     handle = figure;
%     loglog(kp, em_avg + ev_avg, 'LineWidth', 1.25);
%     hold on;
% 
%     loglog(kp, ev_avg);
%     loglog(kp, em_avg);
%     legend('total', 'variance', 'mean');
% 
%   % % % % % % % % % % % % % % % % 
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
% end
% cases = loadcases();
% 
% fsavefigs = true;
% 
% for casename = cases
% 
% 
% % datafolder = fullfile('data',char(casename));
% figfolder = fullfile('ggg');
% if fsavefigs && ~isdir(figfolder); mkdir(figfolder); end
% 
% dataext = load(fullfile('diag11.mat'));
% 
% qhc = coarsenqh('', 0); % load data, second argument is junk
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
% 








% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Ahamd Branch
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


% cases = loadcases();
% 
% fsavefigs = true;
% 
% for casename = cases
% 
% 
% figfolder = fullfile('ggg');
% if fsavefigs && ~isdir(figfolder); mkdir(figfolder); end
% 
% dataext = load(fullfile('diag11.mat'));
% 
% qhc = coarsenqh('', 0); % load data, second argument is junk
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