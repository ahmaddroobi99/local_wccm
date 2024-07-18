fprintf('initializing path directories...\n');

filename = mfilename('fullpath');
idx = max(strfind(filename,filesep));
dirname = filename(1:idx);

addpath(dirname);
addpath(genpath(fullfile(dirname,'solvers')));
addpath(genpath(fullfile(dirname,'drivers')));
addpath(genpath(fullfile(dirname,'utils')));
addpath(genpath(fullfile(dirname,'analysis')));

clear idx filename dirname
