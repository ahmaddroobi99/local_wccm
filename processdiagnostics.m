clear;  
fsavefigs = 1;

cases = loadcases();
casename = cases(3);
datarootfolder = 'data';
datafolder = fullfile(datarootfolder,casename);
parts = strsplit(datafolder, filesep);
name = parts{2};
figfolder = fullfile('processdiagnostics_fig');
if fsavefigs && ~isdir(figfolder); mkdir(figfolder); end

filterdiagfname = fullfile('filterdiagnostics.txt');


T = readtable(filterdiagfname);
T = unique(T);

Ltracersa = unique(T.Ltracers);
sig0a = unique(T.sigma_o);

handle = figure(1); clf;

for j = 1:length(Ltracersa)
    
    for i = 1:length(sig0a)
        
        sig0 = sig0a(i);
        
        matchedinds = (T.sigma_o == sig0) & (Ltracersa(j) == T.Ltracers);
        
        Tn = unique(T(matchedinds,:));
        
        subplot(1,3,1);
        l1(i) = plot(Tn.N,Tn.rms_ave, '.-');
        hold on;
        
        xlabel('(2N+1)');
        ylabel('RMS')
        
        subplot(1,3,2);
        l2(i) = plot(Tn.N,Tn.nrms_ave, '.-');
        hold on;
        
        xlabel('(2N+1)');
        ylabel('NRMS')
        
        
        subplot(1,3,3);
        l3(i) = plot(Tn.N,Tn.xcor_ave, '.-');
        hold on;
        
        xlabel('(2N+1)');
        ylabel('XCorr');
        
        
    end
    legend(l3,arrayfun(@(x) sprintf('\\sigma_o = %.2g',x),sig0a,'UniformOutput',false))
    suptitle(sprintf('tracers L = %i', Ltracersa(j)));
    
    drawnow;
    set(gcf,'position',[1087 552 816 300]);
    export_fig(fullfile(figfolder, sprintf('rmsxcorr_asN-L%i', Ltracersa(j))), '-q101','-p.01','-m1','-pdf');
    savefig(handle,fullfile(figfolder, sprintf('rmsxcorr_asN-L%i', Ltracersa(j))),'compact')
    
end


T = readtable(filterdiagfname);
T = unique(T);

Na = unique(T.N);
sig0a = unique(T.sigma_o);

handle = figure(1); clf;

for j = 1:length(Na)
    
    for i = 1:length(sig0a)
        
        sig0 = sig0a(i);
        
        matchedinds = (T.sigma_o == sig0) & (Na(j) == T.N);
        
        Tn = unique(T(matchedinds,:));
        
        subplot(1,3,1);
        l1(i) = plot(Tn.Ltracers,Tn.rms_ave, '.-');
        hold on;
        
        xlabel('L');
        ylabel('RMS')
        
        subplot(1,3,2);
        l2(i) = plot(Tn.Ltracers,Tn.nrms_ave, '.-');
        hold on;
        
        xlabel('L');
        ylabel('NRMS')
        
        
        subplot(1,3,3);
        l3(i) = plot(Tn.Ltracers,Tn.xcor_ave, '.-');
        hold on;
        
        xlabel('L');
        ylabel('XCorr');
        
        
    end
    legend(l3,arrayfun(@(x) sprintf('\\sigma_o = %.2g',x),sig0a,'UniformOutput',false))
    suptitle(sprintf(' (2N+1)^2 = %i^2', Na(j)));    
    
    drawnow;
    set(gcf,'position',[1087 552 816 300]);
    export_fig(fullfile(figfolder, sprintf('rmsxcorr_asL-N%i', Na(j))), '-q101','-p.01','-m1','-pdf');
    savefig(handle,fullfile(figfolder, sprintf('rmsxcorr_asL-N%i', Na(j))),'compact')
    
end



