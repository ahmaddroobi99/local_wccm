% clear;

% datafolderroot = 'data';
%
% set(1).datafolder = fullfile(datafolderroot,'baroRK4_N128_beta0_d0.005_nu2.5e-28_forcing-sin_tau0-0.04_dt0.02_Nt30000_IC-supplied_set-02-27-18-212412');
% set(1).sigmax = 0.15;
%
% set(2).datafolder = fullfile(datafolderroot,'baroRK4_N128_beta0_d0.005_nu6.3e-25_forcing-constant_tau0-0_kf30_epsf1_sigf1_dt0.001_Nt30000_IC-supplied_set-02-27-18-212412');
% set(2).sigmax = 25;
%
% fsavedata = true;
% Ltracers_  = [16 26 72 144];
%
% kk = 1;
% for nn = 1:length(set)
%     for i = 1:length(Ltracers_)
%         kk = 1 + kk;
%         F(kk) = parfeval(@generate_baro_obs, 0, set(nn).datafolder, Ltracers_(i), set(nn).sigmax, 'savedata', true);
%     end
% end
%
%
% while ~all(arrayfun(@(x) strcmp(x, 'finished'), {F.State}))
%     pause(10)
%     F.Diary
% end

clear;

% datafolderroot = 'data';
% cases = loadcases();

% set = struct();
% for i = 1:length(cases)
%     set(i).datafolder = fullfile(datafolderroot,char(cases(i)));
% end

fsavedata = true;
Ltracers_  = [10 20 40 80];
sigmax_ = [0.005, 0.01, 0.1, 0.2];
% sigmax_=[0];

for i = 1:length(Ltracers_)
    fprintf('running cases with L_tracers = %i\n',Ltracers_(i));
    
    for j = 1:length(sigmax_)
        fprintf('\trunning cases with sigma0 = %.2g\n', sigmax_(j));
        % for k = 1:length(set)
            F = parfeval(@generate_baro_obs, 1, '', Ltracers_(i), sigmax_(j), 'savedata', true);
        % end
%         while ~all(arrayfun(@(x) strcmp(x, 'finished'), {F.State}))
%             pause(1)
%             disp(F(1).Diary);
%         end
        fetchOutputs(F);
%         afterEach(F, fprintf('\t\tdone a case'), 0);
    end
end



