function UVnovo()


% Add directories associated with UVnovo to the Matlab search path.
uvnovo_root = [fileparts(mfilename('fullpath')), filesep];
append_dirs = {'','core_ms','core_ms\fileIO','dependencies', ...
			   'UVnovo','UVnovo\deprecate','UVnovo\models','utils'};
addpath( strjoin( strcat(uvnovo_root, append_dirs), ';') )

%%

% @TODO port an old ui over here.

%%	Prompt for execution mode
%	(1) Data initilization: Import data and create training and test sets.
%	(2) UVnovo training: construct random forests from training data.
%	(3) De novo sequencing: perform de novo sequencing on unknown spectra.
%	(4) Benchmarking: compare de novo results to PSM knowledge.

% UVnovo_init_data(spectra_file, psm_file)
UVnovo_init_data()







