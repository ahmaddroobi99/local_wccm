function cstr = casename(N,Ltracers,sigmax,Tend,spectrum_type,damping_type,dt)

cstr = sprintf('N%i_L%i_sigmao%.2g_Te%0.0f_sp-%s_da-%s_dt%.3g%', N, Ltracers, sigmax, Tend, spectrum_type, damping_type, dt);

end