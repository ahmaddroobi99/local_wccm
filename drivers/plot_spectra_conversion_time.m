% cases = 's';

% fsavefigs = true;
% set(gcf,'Color','w');

% for casename = cases
    
    
% datafolder = fullfile('data',char(casename));
% figfolder = fullfile('work_on_spectra_converion', char(casename));
% if fsavefigs && ~isdir(figfolder); mkdir(figfolder); end

% dataext = load(fullfile('diag11.mat'));

% qhc = coarsenqh(datafolder, 0.02); % load data, second argument is junk

% Ntruth = sqrt(size(qhc,1));
% k = [0:Ntruth/2 -Ntruth/2+1:-1]';
% [KX, KY] = meshgrid(k,k);
% Knorm = hypot(KX,KY);
% Knorminv = 1./Knorm; Knorminv(1,1) = 0;

% kp = Knorm(1,1:Ntruth/2+1);


% psi_hat_qhc = -Knorminv(:).*qhc/Ntruth^2;
% t_columns = size(psi_hat_qhc,2);   


% % display(t_columns);
% t = 2*10^5;
% dt = 0.02;
% % t/dt;
% tslice= 1 ;

% num_slices = 1;  

% all_em_ev = zeros(length(kp), 2, tslice);

% for slice_idx = 1:tslice

%     psi_hat_qhc_slice = psi_hat_qhc;
%     qhcm = abs(mean(psi_hat_qhc_slice,2)).^2;
%     qhcv = var(psi_hat_qhc_slice,0,2);

%     qhcm = reshape(qhcm,Ntruth,Ntruth);
%     qhcv = reshape(qhcv,Ntruth,Ntruth);

%     ev = zeros(Ntruth/2+1,1);
%     em = zeros(Ntruth/2+1,1);
%     for jj = 1:Ntruth
%         for ii = 1:Ntruth
%             k = sqrt(KX(ii,jj).^2 + KY(ii,jj).^2);
%             if ceil(k) <= Ntruth/2
%                 r = k-floor(k);
%                 ev(floor(k)+1) = ev(floor(k)+1) + (1-r)*abs(qhcv(ii,jj));
%                 ev(ceil(k)+1) = ev(ceil(k)+1) + r*abs(qhcv(ii,jj));

%                 em(floor(k)+1) = em(floor(k)+1) + (1-r)*abs(qhcm(ii,jj));
%                 em(ceil(k)+1) = em(ceil(k)+1) + r*abs(qhcm(ii,jj));
%             end
%         end
%     end
%     ev = 0.5*ev;
%     em = 0.5*em;

%     all_em_ev(:,:,slice_idx) = [em, ev];
% end
% handle = figure;

% hold on;

% for slice_idx = 1:num_slices
%     em = all_em_ev(:,1,slice_idx);
%     ev = all_em_ev(:,2,slice_idx);
%     loglog(kp, em + ev, 'LineWidth', 1.25);
% end


% if strcmp(dataext.params.forcingtype,'constant')
%     kkp = (1:35)';

%     alpha = -5/3;
% else
%         kkp = (1:45)';

%     alpha = -3;
% end

% [~,~,ikp] = intersect(kkp,kp);
% amin = fminsearch(@(a) mean(log(abs(a*kkp.^alpha - ev(ikp)))), 1);

% loglog(kkp,amin*kkp.^alpha,'k--');
% text(2,amin*kkp(1).^alpha,strcat('k^{',strtrim(rats(alpha)),'}'))

% xlabel('wavenumber');
% ylabel('energy');
% legend('Slice 1', 'Slice 2', 'Slice 3', 'Slice 4', 'Slice 5');




% if fsavefigs
%     resizefonts(handle);
%     export_fig(fullfile(figfolder, 'Hi testing'), '-q101','-p.01','-m1','-pdf');
%     % savefig(handle, fullfile(figfolder, 'all_slices_spectrum.fig'));
% end
% end
% % 555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555

% % psi_hat_qhc =psi_hat_qhc(:,1:tslice);
% % qhcm = abs(mean(psi_hat_qhc,2)).^2;
% % qhcv = var(psi_hat_qhc,0,2);

% % qhcm = reshape(qhcm,Ntruth,Ntruth);
% % qhcv = reshape(qhcv,Ntruth,Ntruth);

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
            
% %             em(floor(k)+1) = em(floor(k)+1) + (1-r)*abs(qhcm(ii,jj));
% %             em(ceil(k)+1) = em(ceil(k)+1) + r*abs(qhcm(ii,jj));
% %         end
% %     end
% % end
% % ev = 0.5*ev;
% % em = 0.5*em;


% % %%
% % handle = figure;
% % loglog(kp,em+ev,'LineWidth',1.25);
% % hold on;

% % loglog(kp,ev); 
% % loglog(kp,em);
% % legend('total','variance','mean');


% % if strcmp(dataext.params.forcingtype,'constant')
% %     kkp = (1:35)';

% %     alpha = -5/3;
% % else
% %         kkp = (1:45)';

% %     alpha = -3;
% % end
% % [~,~,ikp] = intersect(kkp,kp);
% % amin = fminsearch(@(a) mean(log(abs(a*kkp.^alpha - ev(ikp)))), 1);

% % loglog(kkp,amin*kkp.^alpha,'k--');
% % text(2,amin*kkp(1).^alpha,strcat('k^{',strtrim(rats(alpha)),'}'))



% % xlabel('wavenumber');
% % ylabel('energy')

% % drawnow; 
% % if fsavefigs
% %     resizefonts(handle);
% %     export_fig(fullfile(figfolder, sprintf('spectrum')), '-q101','-p.01','-m1','-pdf');
% %     savefig(handle, fullfile(figfolder, sprintf('spectrum.fig')))
% % end



% % end




% 10/3/2024



cases = 'm';

fsavefigs = true;

for casename = cases
    
    
% datafolder = fullfile('data',char(casename));
figfolder = fullfile('ssss', char(casename));
if fsavefigs && ~isdir(figfolder); mkdir(figfolder); end

dataext = load(fullfile('diag11.mat'));

qhc = coarsenqh('', 0.02); % load data, second argument is junk

Ntruth = sqrt(size(qhc,1));
k = [0:Ntruth/2 -Ntruth/2+1:-1]';
[KX, KY] = meshgrid(k,k);
Knorm = hypot(KX,KY);
Knorminv = 1./Knorm; Knorminv(1,1) = 0;

kp = Knorm(1,1:Ntruth/2+1);


psi_hat_qhc = -Knorminv(:).*qhc/Ntruth^2;
qhcm = abs(mean(psi_hat_qhc,2)).^2;
qhcv = var(psi_hat_qhc,0,2);

qhcm = reshape(qhcm,Ntruth,Ntruth);
qhcv = reshape(qhcv,Ntruth,Ntruth);

ev = zeros(Ntruth/2+1,1);
em = zeros(Ntruth/2+1,1);
for jj = 1:Ntruth
    for ii = 1:Ntruth
        k = sqrt(KX(ii,jj).^2 + KY(ii,jj).^2);
        if ceil(k) <= Ntruth/2
            r = k-floor(k);
            % plus 1 because of matlab indexing
            ev(floor(k)+1) = ev(floor(k)+1) + (1-r)*abs(qhcv(ii,jj));
            ev(ceil(k)+1) = ev(ceil(k)+1) + r*abs(qhcv(ii,jj));
            
            em(floor(k)+1) = em(floor(k)+1) + (1-r)*abs(qhcm(ii,jj));
            em(ceil(k)+1) = em(ceil(k)+1) + r*abs(qhcm(ii,jj));
        end
    end
end
ev = 0.5*ev;
em = 0.5*em;


%%
handle = figure;
loglog(kp,em+ev,'LineWidth',1.25);
hold on;

loglog(kp,ev); 
loglog(kp,em);
legend('total','variance','mean');


if strcmp(dataext.params.forcingtype,'constant')
    kkp = (1:35)';

    alpha = -5/3;
else
        kkp = (1:45)';

    alpha = -3;
end
[~,~,ikp] = intersect(kkp,kp);
amin = fminsearch(@(a) mean(log(abs(a*kkp.^alpha - ev(ikp)))), 1);

loglog(kkp,amin*kkp.^alpha,'k--');
text(2,amin*kkp(1).^alpha,strcat('k^{',strtrim(rats(alpha)),'}'))



xlabel('wavenumber');
ylabel('energy')

drawnow; 
if fsavefigs
    resizefonts(handle);
    export_fig(fullfile(figfolder, sprintf('spectrum')), '-q101','-p.01','-m1','-pdf');
    savefig(handle, fullfile(figfolder, sprintf('spectrum.fig')))
end



end