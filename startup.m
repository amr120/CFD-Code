% Set default figure position
set(0, 'DefaultFigurePosition', [2195 403 708 420])

% Load all folders and subfolders on the default path
paths = genpath(pwd); paths = regexp(paths,':','split');
paths = paths(contains(paths,{'Backup','Figures','Export','evtol'}) == 0);
path_list = [];
for n = 1:length(paths); path_list = [path_list paths{n} ':']; end;
addpath(path_list);
