% coarsen raw data set


function [params] = train_msm(datafolder,dtsave)

savefilename = fullfile(datafolder,sprintf('trainedmsm_dt%.2g.mat',dtsave));

searchfilename = dir(fullfile(datafolder, sprintf('trainedmsm_dt%.2g.mat',dtsave))); % arbitrarily load regardless of dtsave size
if ~isempty(searchfilename)
    data = load(fullfile(searchfilename.folder,searchfilename.name));
    params = data.params;
    return
end


fnameqh = fullfile('qh11.bin');
fnamet = fullfile('t11.bin');

fidqh = fopen(fnameqh,'r');
fidt = fopen(fnamet,'r');

t = fread(fidt,'double');
n_t = length(t);
dt = t(2) - t(1);
Tend = t(end);

dtsave = round(dtsave/dt)*dt; % round to a factor of dt
nsavef = dtsave/dt;
tsave = t(1:nsavef:length(t));
n_tsave = length(tsave);

dataext = load(fullfile('diag11.mat'));
Ntruth =  dataext.params.N;
Ntotal = Ntruth*Ntruth;


qhc = zeros(Ntotal, n_tsave);
ii = 0;
for i = 1:n_t
    qq = fread(fidqh,2*Ntotal,'double');
    qq = complex(qq(1:Ntotal), qq(Ntotal+1:2*Ntotal));
    if ~mod(i-1,nsavef)
        ii = ii + 1;
        qhc(:,ii) = qq(:);
    end
end
fclose(fidqh); fclose(fidt);


N = sqrt(size(qhc,1));
k = [0:Ntruth/2 -Ntruth/2+1:-1]';
[KX, KY] = meshgrid(k,k);
Knorm = hypot(KX,KY);
Knorminv = 1./Knorm; Knorminv(1,1) = 0;
Ut = -Knorminv(:).*qhc/Ntruth^2;
Ek = var(Ut,0,2);

d = zeros(Ntotal,1);
omega = zeros(Ntotal,1);
sigma = zeros(Ntotal,1);

for i = 1:Ntruth^2
    if ~mod(i,100)
    fprintf('%.2f\n',i/Ntruth^2)
    end
    [Rmm,lags] = xcorr(Ut(i,1:end).'- mean(Ut(i,1:end),2), 'coeff');
    Rmm = Rmm(lags>0);
    lags = lags(lags>0);
    
    
    
    steady = colfilt(Rmm.', [1, round(25/dt)], 'sliding', @(x) all(abs(real(x)) <= 0.1));
    indc = find(diff(steady) == 1);
    
    truncatedR = Rmm(1:indc);
%     clf;
%     plot(lags(1:indc)*dt, ( real(truncatedR)));
%     drawnow;
    
    rint = sum(truncatedR)*dt;
    
    Tk = real(rint);
    if Tk < 0
        Tk = 0.01; % small value
    end
    thetak = -imag(rint);
    
    d(i) = Tk/(Tk^2+thetak^2);
    omega(i) = thetak/(Tk^2+thetak^2);
    sigma(i) = sqrt(2*Ek(i)*d(i));
    
end

params.d = d;
params.omega = omega;
params.sigma = sigma;
save(savefilename,'params','dt','Tend');

end