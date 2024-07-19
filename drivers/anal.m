function T = analysis(pin, opts, m, qhc, Xt)
if any(any(isnan(m)))
    warning('unstable, nan detected');
    m(isnan(m)) = 0;
    T = table; 
end
casestr = casename(pin.N, pin.Ltracers, pin.sigmax, pin.Tend, pin.spectrum_type, pin.damping_type, pin.dtsave);

N = pin.computed.N;
L = 2*pi;
x = L/N*(0:N-1)';               % x coordinate
y = L/N*(0:N-1)';               % y coordinate

dtsave = pin.computed.dtsave;
Ntotal = pin.computed.Ntotal;
tsave = pin.computed.tsave;
n_tsave = length(tsave);
plot_time = opts.plot_times;

sigma = pin.computed.sigma;
gamma = pin.computed.gamma;
Sigma = pin.computed.Sigma;
Gamma = pin.computed.Gamma;

E0 = pin.E0;
Ltracers = pin.Ltracers;
fsavefigs = opts.fsavefigs;

if ~isfield(pin,'datafolder')
    pin.datafolder = '';
end

datafolder = pin.datafolder;
figfolder = fullfile('new_wccm', casestr);
if fsavefigs && ~isdir(figfolder); mkdir(figfolder); end

ptidx = floor(plot_time./dtsave + 1);
tdecor = 1/min(nonzeros(Gamma)); % decorrelation time of the slowest mode
time_burnin = 1*tdecor; %  discard 1 decor unit of time
idxburnin = floor(time_burnin./dtsave+1);
if idxburnin > n_tsave
    warning('burnin period longer then simulation end time...increase Tend...taking all time units');
    idxburnin = 1;
end
axcolors = get(gca,'Defaultaxescolororder');
Ntruth = sqrt(size(qhc,1));
k = [0:Ntruth/2 -Ntruth/2+1:-1]';
[KX, KY] = meshgrid(k,k);
Knorm = hypot(KX,KY);
Knorminv = 1./Knorm;
Knorminv(1,1) = 0;
Ut = -Knorminv(:).*qhc/Ntruth^2;
Ut(junkmodes2d(Ntruth)) = 0;
Utm = abs(mean(Ut,2)).^2;
Utv = var(Ut,0,2);

Ev_truth_2d = reshape(Utm,Ntruth,Ntruth);
Em_truth_2d = reshape(Utv,Ntruth,Ntruth);
[Ev_truth,kv_truth] = spectra_omnidirectional_sum(Em_truth_2d);
[Em_truth,km_truth] = spectra_omnidirectional_sum(Ev_truth_2d);
E_truth = Em_truth + Ev_truth;

tendfactor = round(tsave(end)/tdecor);
tdecorfactor = unique([10 25 50 100 150  tendfactor],'stable');

if ~all(tdecorfactor<=tendfactor)
    warning(sprintf('Averaging time is %.2g decor factors, you may want to increase simulation time',tendfactor));
    tdecorfactor = tdecorfactor(tdecorfactor<=tendfactor);
end

taveidx = floor(tdecor*tdecorfactor./dtsave + 1);
for jj = 1:length(taveidx)

   indtime = idxburnin:taveidx(jj);
   indsubsamp = subsample2d(Ntruth,N);
   Ev_analytical_full = (sigma.*sigma)./(2*gamma);
   Ev_analytical_full(junkmodes2d(N)) = 0;
   [Ev_analytical, kv_analytical,count] = spectra_omnidirectional_sum(Ev_analytical_full); % mean zero by theory
   Ev_recovery_2d = reshape(var(m(:,indtime), 0, 2), N, N);
   Em_recovery_2d = reshape(abs(mean(m(:,indtime), 2)).^2, N, N);
   [Ev_recovery, kv_recovery] = spectra_omnidirectional_sum(Ev_recovery_2d);
   [Em_recovery, km_recovery] = spectra_omnidirectional_sum(Em_recovery_2d);
   E_recovery = Em_recovery + Ev_recovery;
    %-------------------------------------------------------------------------%
    % plot averaged spectral recovery by wavenumber
    %-------------------------------------------------------------------------%
    handle = figure(1); clf;
    set(handle,'Color','w');
    subplot_tight(1,2,1,[0.11,0.08]);
    loglog(kv_truth, E_truth , '-','color',axcolors(1,:)); hold on;
    loglog(kv_recovery, E_recovery,'-','markersize',12,'color',axcolors(2,:));
    xlabel('wavenumber k');
    ylabel('E_k')
    lg = legend('truth kv', 'truth E ', 'location','southwest');
    lg.FontSize = 8;
    xlim([1 max(kv_truth)])
    ylim([10^-11 5*max(E_truth)]);
    title('total energy') 
    hold on
    
    ax = subplot_tight(1,2,2,[0.11,0.08]);
    ll1=loglog(kv_truth, Ev_truth , '--', 'color',axcolors(1,:)); hold on;
    ll2=loglog(kv_truth, Em_truth , '-', 'color',axcolors(1,:)); hold on;

    ll3=loglog(kv_recovery, Ev_recovery,'--','color',axcolors(2,:));
    ll4=loglog(kv_recovery, Em_recovery,'-','color',axcolors(2,:));
    lg = legend([ll1 ll3, ll2  ll4],{'truth variance','recovery variance',  'truth mean', 'recovery mean'}, 'location','southwest'); 
    lg.FontSize = 8;
    xlabel('wavenumber k');
    xlim([1 max(kv_truth)])
    ylim([10^-11 5*max(E_truth)]);
    title('mean  and variance energy') 

    annotation('textbox', [0 0.9 1 0.1], ...
        'String', 'spectral recovery', ...
        'EdgeColor', 'none', ...
        'HorizontalAlignment', 'center',...
        'FontSize',12)
    set(handle, 'Position', [459 413 674 372]);

    drawnow;
    if fsavefigs
        resizefonts(handle);
        export_fig(fullfile(figfolder, sprintf('spec-decorf%i',tdecorfactor(jj))));
        savefig(handle, fullfile(figfolder, sprintf('spec-decorf%i.fig',tdecorfactor(jj))),'compact')
    end
    
    %-------------------------------------------------------------------------%
    % plot averaged rmse and pattern correlation recovery by wavenumber
    %-------------------------------------------------------------------------%
   vf = m(:,indtime);
   vt = Ut(indsubsamp,indtime);
   rms      = sqrt(sum((vf - vt).* conj(vf - vt), 1));

   rms_mode = sqrt((vf - vt).*conj(vf - vt));
   nrms      = rms./sqrt(sum(vt.*conj(vt), 1));
   nrms_mode = rms./sqrt(vt.*conj(vt));
   xcor      = real( sum(vf.*conj(vt), 1) ./ ( sqrt(sum(vf.*conj(vf), 1)) .* sqrt(sum(vt.*conj(vt), 1))) );
   xcor_mode = real( vf.*conj(vt) ./ ( sqrt(vf.*conj(vf)) .* sqrt(vt.*conj(vt))) );
   
   rms_ave  = mean(rms);
   nrms_ave = mean(nrms);
   rms_mode_ave  = mean(rms_mode,2);
   nrms_mode_ave = mean(nrms_mode,2);
   [rms_mode_ave_omni, knormf]  = omnidirectional_sum(rms_mode_ave);
   [nrms_mode_ave_omni, knormf] = omnidirectional_sum(nrms_mode_ave);

   xcor_ave = mean(xcor(2:end));
   xcor_mode_ave = mean(xcor_mode,2);
   [xcor_mode_ave_omni, knormf] = omnidirectional_sum(xcor_mode_ave);

   
   
    handle = figure(3); clf;
    set(handle,'Color','w');

    yyaxis left
    plot(knormf, rms_mode_ave_omni );
    hold on;

    xlim([1 max(knormf)])
    plot([1 max(knormf)], pin.sigmax*[1 1])
    hold off;

    ylabel('RMSE');
    ylimc = ylim;
    ylim([0 ylimc(2)*1.1]);
    text(1.5,pin.sigmax*1.05,'obs noise \sigma^o','color',axcolors(1,:)) 
    
    yyaxis right
    plot(knormf, nrms_mode_ave_omni );
    ylabel('NRMSE');
    xlim([1 max(knormf)])
    ylimc = ylim;
    ylim([0 ylimc(2)*1.1]);
    xlabel('wavenumber k');
    set(handle, 'Position', [799 137 399 372]);
    title('Average RMS error')

    drawnow;
    if fsavefigs
        resizefonts(handle);
        export_fig(fullfile(figfolder, sprintf('rmsmode-decorf%i',tdecorfactor(jj))), '-q101','-p.01','-m1','-pdf');
        savefig(handle, fullfile(figfolder, sprintf('rmsmode-decorf%i.fig',tdecorfactor(jj))),'compact')
    end

    handle = figure(3);clf;
    plot(knormf, xcor_mode_ave_omni );
    xlim([1 max(knormf)])
    ylabel('XCorr');
	xlabel('wavenumber k');
    ylim([0 1]);
    
    set(handle, 'Position', [799 137 399 372]);
    title('Average XCorr')
    drawnow;

    if fsavefigs
        resizefonts(handle);
        export_fig(fullfile(figfolder, sprintf('xcorrmode-decorf%i',tdecorfactor(jj))), '-q101','-p.01','-m1','-pdf');
        savefig(handle, fullfile(figfolder, sprintf('xcorrmode-decorf%i.fig',tdecorfactor(jj))),'compact')
    end
        

    
    %-------------------------------------------------------------------------%
    % plot spectral recovery surface
    %-------------------------------------------------------------------------%
    
    handle = figure(4);
    set(handle,'Color','w');
    clf; 
    E_truth_2d_log = fftshift(log10(Ev_truth_2d + Em_truth_2d));
    E_truth_2d_log = E_truth_2d_log(2:end,2:end);
    KXp = fftshift(KX);
    KYp = fftshift(KY);
    KXp = KXp(2:end,2:end);
    KYp = KYp(2:end,2:end);    
    E_truth_2d_log(isinf(E_truth_2d_log)) = 0;
    
    
    E_recovery_2d_log = fftshift(log10(Ev_recovery_2d + Em_recovery_2d));
    E_recovery_2d_log = E_recovery_2d_log(2:end,2:end);
    kxxp = fftshift(pin.computed.kxx);
    kyyp = fftshift(pin.computed.kyy);
    kxxp = kxxp(2:end,2:end);
    kyyp = kyyp(2:end,2:end);
    set(gcf,'Position', [268 723 568 271]);
    ax = subplot(1,2,1);
    surf(KXp, KYp, E_truth_2d_log);
    xlabel('k_x');
    ylabel('k_y');
    zlabel('log10(E_k)')
    shading interp;  box on; axis tight;
    title('truth');
    
    s3 = subplot(1,2,2);
    surf(kxxp, kyyp, E_recovery_2d_log)
    xlabel('k_x');
    ylabel('k_y');
    zlabel('log10(E_k)')
    shading interp;  box on; axis tight;
    title('recovery');
    s3Pos = get(s3,'position');
    ch = colorbar('location','eastoutside', 'FontSize', 8);
    set(ch, 'Position', [0.94 0.1100 0.014 0.7150]);
    set(s3,'position',s3Pos);
    suptitle('total energy spectrum')
    
    drawnow;
    if fsavefigs
        resizefonts(handle);
        export_fig(fullfile(figfolder, sprintf('spec2d-decorf%i',tdecorfactor(jj))), '-q101','-p.01','-m2','-png');
        savefig(handle, fullfile(figfolder, sprintf('spec2d-decorf%i.fig',tdecorfactor(jj))),'compact')
    end
   end 
    handle = figure(3); 
    set(handle,'Color','w');
    clf;
    subplot(2,1,1);
    plot(tsave(idxburnin:end),rms);
    hold on;
    a1=plot([tsave(idxburnin) tsave(end)],rms_ave*[1 1],'color','k','linewidth',1.25);
    a2=plot([tsave(idxburnin) tsave(end)],pin.sigmax*[1 1],'--','color','k','linewidth',1.25);
    xlim([tsave(idxburnin) tsave(end)])
    legend([a1,a2],'average RMSE','obs noise \sigma^o');
    ylabel('RMSE');
    title('RMS and XCorr');
    ylimc = ylim;
    ylim([0 ylimc(2)*1.1]);
    
    subplot(2,1,2);
    plot(tsave(idxburnin:end),xcor);
    hold on;
    a1=plot([tsave(idxburnin) tsave(end)],xcor_ave*[1 1],'color','k','linewidth',1.25);
    xlim([tsave(indtime(1)) tsave(end)])
    legend(a1,'average XCorr');
    ylabel('XCorr');
    xlabel('t');
    ylimc = ylim;
    ylim([0 1]);
    set(handle, 'Position', [997 426 620 429]);
    if fsavefigs
        resizefonts(handle);
        export_fig(fullfile(figfolder, sprintf('rmsxcorr-decorf%i',tdecorfactor(jj))), '-q101','-p.01','-m1','-pdf');
        savefig(handle, fullfile(figfolder, sprintf('rmsxcorr-decorf%i.fig',tdecorfactor(jj))),'compact')
    end
    %-------------------------------------------------------------------------%
    % plot true and recovered velocity field and streamlines together
    %-------------------------------------------------------------------------%
    handle = figure(5); clf
    set(gcf,'Position',[440 486 724 346]);
    set(handle,'Color','w');
    XX = zeros(2, n_tsave);
    for i = 1:length(plot_time)
        clf;
        subplot_tight(1,2,1, [0.1, 0.07]);
        [xxq, yyq, Psiqt] = streamfunction(reshape(Ut(:,ptidx(i)), Ntruth, Ntruth),2);
        contourf(xxq, yyq, real(Psiqt), 25); view(2); shading interp;  box on; axis tight;
        clims = caxis;
        colorbar('location','southoutside');
        set(gca,'XTick',[]);
        set(gca,'YTick',[]);
        
        xlim([0 x(end)]);
        ylim([0 x(end)]);
        set(gca,'Layer','top');
        title('streamfunction truth');
        hold on;
        
        for n = 1:Ltracers
            XX(1,:) = mod(Xt((n-1)*2 + 1,:), x(end));
            XX(2,:) = mod(Xt((n-1)*2 + 2,:), x(end));
            plot(XX(1,ptidx(i)), XX(2,ptidx(i)), '.',  'MarkerSize', 12)
        end
    
        s2 = subplot_tight(1,2,2, [0.1, 0.07]);
        [xxq, yyq, Psiqr] = streamfunction(reshape(m(:,ptidx(i)), N, N),4);
        contourf(xxq, yyq, real(Psiqr), 25); view(2); shading interp;  box on; axis tight;
        set(gca,'XTick',[]);
        set(gca,'YTick',[]);
        caxis(clims);
        hold on;
        s2Pos = get(s2,'position');
        ch = colorbar('location','eastoutside', 'FontSize', 8);
        
        set(s2,'position',s2Pos);
        chpos = get(ch,'Position');
        set(ch,'Position',[chpos(1) chpos(2) 0.02 0.7746]);
        
        xlim([0 x(end)]);
        ylim([0 x(end)]);
        set(gca,'Layer','top');
        title('streamfunction recovery');
        for n = 1:Ltracers
            XX(1,:) = mod(Xt((n-1)*2 + 1,:), x(end));
            XX(2,:) = mod(Xt((n-1)*2 + 2,:), x(end));
            plot(XX(1,ptidx(i)), XX(2,ptidx(i)), '.',  'MarkerSize', 12)
        end
    
        suptitle(sprintf('T = %i', round(plot_time(i))));
        drawnow limitrate;
    end
    sigma_o = pin.sigmax;
    tmax = tsave(end);
    T = table(N,rms_ave,nrms_ave,xcor_ave,sigma_o,dtsave,tmax,Ltracers);