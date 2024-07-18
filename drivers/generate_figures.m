%-------------------------------------------------------------------------%
% recursively find .fig files and generate nice pdfs from them
%-------------------------------------------------------------------------%

function generate_figures(folder,recursive)

if nargin < 1
    folder = 'figs'; % by default assume figs is the default folder anme
    if ~isfolder(folder)
        folder = '.'; % if figs is not the default folder assume the current working directory is the default folder
    end
end

if nargin < 2
    recursive = true;
end

if ~isfolder(folder)
    error('folder not found')
end

% Get all .fig files in the folder
files = dir(fullfile(folder, '*.fig'));
files = fullfile(folder, {files.name});

for i = 1:numel(files)
    fprintf('processing %s\n',files{i});

    % Get the filename
    [~, fname] = fileparts(files{i});
    
    % Open and display the .fig file
    hfig = openfig(files{i});
    
    % RESIZE FONTS
    resizefonts(hfig)
    
    % Save as file with the same name as the .fig file
    try
    export_fig(hfig, fullfile(folder, fname), '-q101','-p.01','-m1','-png');
    catch
        warning('failed')
    end
    
    % Close the figure again
    close(hfig)
end


dirs = dir(folder);
dirs(~[dirs.isdir]) = []; %remove non-directories
dirs(ismember( {dirs.name}, {'.', '..'})) = []; %remove current and parent directory.

if recursive
    for i = 1:numel(dirs)
        generate_figures(fullfile(folder, dirs(i).name), recursive)
    end
end