clear;

fsavedata = true;
Ltracers_ = [10];
sigmax_ = [0.005];

for i = 1:length(Ltracers_)
    fprintf('running cases with L_tracers = %i\n', Ltracers_(i));
    
    for j = 1:length(sigmax_)
        fprintf('\trunning cases with sigma0 = %.2g\n', sigmax_(j));
        
        % Call the function with correct arguments
        generate_eulerian_obs('', Ltracers_(i), sigmax_(j), 'savedata', fsavedata);
        
        fprintf('\t\tdone a case\n');
    end
end