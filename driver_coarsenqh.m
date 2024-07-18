
cases = loadcases();



for casename = cases
    
    datafolder = fullfile('data',char(casename));
    dataext = load(fullfile(datafolder, 'diagnostics.mat'));
    fnamet = fullfile(datafolder, 't.bin'); fidt = fopen(fnamet,'r');
    t = fread(fidt,'double'); fclose(fidt);
    dt = t(2) - t(1);
    
    qhc = coarsenqh(datafolder, dt*10);
end