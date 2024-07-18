% %  find  size of truth and observations
% your_matrix_variable= load('obs_Ltracers10_sigmax0.1');  % Adjust the variable name if it differs
% size(your_matrix_variable.X{1})  % Replace 'your_matrix_variable' with the actual variable name


dataext = load(fullfile('diag11.mat'));
Ntruth =  dataext.params.N;
Ntotal = Ntruth*Ntruth;

fnameqh = fullfile('qh11.bin');
fnamet = fullfile('t11.bin');

fidqh = fopen(fnameqh,'r');
qq = fread(fidqh,Ntotal,'double');
fidt = fopen(fnamet,'r');
t = fread(fidt,'double'); 
n_t = length(t);
dt = t(2) - t(1);
display(size(qq));


qhc = coarsenqh('',0.2);

% display(qq) ;
