function params = user_params()
% User defined parameters.
%	Any specified here replace the values located in ./UVnovo/default_params().


%% Data import and setup of training and validation datasets

%%% Import spectral data.
params.pre.import_spectra = struct( ...
	... primary parent mass field name
	'pmass_field', 'pmass_theorNormmass', ...
	... MS2 peak filtering
	'filterPeaks', struct( ...
		'minInten', 5, ...			% remove peaks below this intensity
		'removeZeroInten', true ...	% remove peaks with zero intensity
		) ...
	);

%%% Create cross-validation datasets.
params.pre.crossVal.nPartitions = 3; % number of CV partitions

