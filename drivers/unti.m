    
    datafolder ='Ahmad_QG\N16_L10_sigmao0.1_Te60_sp-kolmo_da-unif_dt0.2\qgahmad_data.mat';
    load(datafolder, 'pin', 'm', 'Xt');
    timer = tic;
    qhc = coarsenqh(datafolder,pin.dtsave);
    fprintf('coarsenqh time %s\n', time2str(toc(timer)))
    opts.renderer = 'zbuffer';
    opts.plot_times = [0:1:60];
    opts.fsavefigs = true;
    anal(pin, opts, m, qhc, Xt)

