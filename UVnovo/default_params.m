function params = default_params()


%% Data import and setup of training and validation datasets

%%% Import spectral data.
params.pre.import_spectra = struct( ...
	'pmass_field', 'pmass', ... % primary parent mass field name
	'filterPeaks', struct( ...		% MSfilterScans()
		'minInten', 5, ...			% remove peaks below this intensity
		'removeZeroInten', true ...	% remove peaks with zero intensity
		) ...
	);

%%% Create cross-validation datasets.
params.pre.crossVal.nPartitions = 3; % number of CV partitions

