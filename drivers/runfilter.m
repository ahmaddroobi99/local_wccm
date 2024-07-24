    cases = loadcases();
    casename = cases(3);

    datafolder = fullfile('data',casename);
    dataext = load(fullfile('diag11.mat'));
    fnamet = fullfile('t11.bin'); fidt = fopen(fnamet,'r');

    t = fread(fidt,'double'); fclose(fidt);

    dt = t(2) - t(1);
    nt = length(t);
    Ntruth =  dataext.params.N;

    Ltracers = 20;
    sigmax = 0.01;

    obsdatafile = dir(fullfile(sprintf('obs_Ltracers%i_sigmax%.2g.mat', Ltracers, sigmax)));
    if isempty(obsdatafile)
        error('obs data not found');
    end

    obsdata = load(fullfile(obsdatafile.folder,obsdatafile.name));
    % pin.Tend = t(end);

    pin.Tend = 30;
    pin.nt=pin.Tend\dt;
    pin.Ltracers = Ltracers;
    pin.N = 16; % Grid points in each dimension (must be even)
    pin.X = obsdata.X;
    pin.E0 = 0.8;
    pin.damping_type = 'unif';
    pin.d = 0.1;
    pin.mu = 0.05;
    pin.spectrum_type = 'kolmo';
    pin.num_modes_flat = 2;
    pin.sigmax = obsdata.sigmax;
    pin.dt = dt;
    pin.dtsave = dt * 10;
    pin.fsavedata = true;
    pin.datafolder = 'AhmadH_QG';
    pin.beta =  0.05;
    pin.dd = 5E-3;
    pin.epsf = 0;
    pin.kf = 50;  % no effect since epsf = 0
    pin.sigf = 1; % no effect since epsf = 0
    pin.alpf = 0; % no effect since epsf = 0
    pin.tau0 = 0.04;
    pin.tdel = 0;
    pin.forcingtype = 'sin';
    pin.Nt = 2e5;
    pin.nu=1;
    pin.sh=0.2;

    timer = tic;
    display("WHY filter begins");
    [pin, m, Xt] = where(pin);
    display("filter ends");
    fprintf('filter runtime  %s\n', time2str(toc(timer)));
    
    timer = tic;
    qhc = coarsenqh(datafolder,pin.dtsave);
    fprintf('coarsenqh time %s\n', time2str(toc(timer)))

    opts.renderer = 'zbuffer';
    opts.plot_times = [0:1:30];
    opts.fsavefigs = true;
    analysis(pin, opts, m, qhc, Xt);
