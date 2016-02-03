function params = default_params()

%% General parameters
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

% Print info messages? No messages [0], some [1], all [2].
% @TODO get this widely implemented or completely remove.
params.UVnovo.verbosity = 2;

% UVnovo training and de novo speed improves greatly with matlabpool open. This
% requires the Matlab parallel computing toolbox.
params.UVnovo.useParallel = true;
% Argument(s) for matlabpool() if using parallel and no pool is open.
params.UVnovo.matlabpool_varargs = {4};


% Amino acids
% @TODO I'm not happy with how this works. Make it better.
params.AAs = struct(  ...
	'excludeAAs', '', ...  % disallowed amino acids
	'ncderiv', [0, 0], ... % N/C-terminal derivatization fixed mass
	'ptms', ... % PTMs: <2 x n> cell of { symbol<char>, mass<real>; ... }
		{{
			'm', CONSTS.mAA.M + CONSTS.mPTM.Oxidation;
		}} ....
	);


%% Data import and setup of training and validation datasets
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

% Import spectral data. Accessed by import_spectra().
params.pre.import_spectra = struct( ...
	'pmass_field', 'pmass', ... % primary parent mass field name
	'filterPeaks', struct( ...	% MSfilterScans()
		'minInten', 5, ...			% remove peaks below this intensity
		'removeZeroInten', true ...	% remove peaks with zero intensity
		) ...
	);

% Define cross-validation datasets.
params.pre.crossVal.nPartitions = 3; % number of CV partitions


%% UVnovo training. Construct random forests from spectra & psms.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

% Mass window, around a base peak, from which to draw feature variables.
params.train.windowRange_Da = [-50 50];

% Seed for repeatable generation of 'random' true-negative training examples.
% If left empty, UVnovo assigns a random seed, and the training procedure will
% yield slightly different results each time it is run.
params.train.falseObservationSeed = 121212;


% Ensemble training schedule. Iteratively train & select top predictor vars for
% each subsequent model.
% Parameters for each training round:
%	nvars: Select this many predictor vars. First step must be '0' (all vars).
%	ntrees: Number of decision trees in the ensemble.
%	plotOOBError: Calculate & plot out-of-bag error vs. ntrees. It can be slow.
%	keepEns: Retain & return specified ensemble model(s).
%	nPrint: Print training progress after every n trees.
%	minLeaf: Minimum number of observations per tree leaf.
%	oobVarImp: <'on'|'off'> Out-of-bag estimate of var importance. SLOW!
params.train.ensemble.iters = cell2struct({
	0,      150,     1,             0,        4,       1,       'off'; % round 1
	200,    150,     1,             0,        4,       1,       'off'; % round 2
	100,    150,     1,             0,        4,       1,       'off'; % round 3
	60,     250,     1,             0,        4,       1,       'on' ; % round 4
	30,     400,     1,             1,        4,       1,       'on' ; % final
	}, ...
  {'nvars','ntrees','plotOOBError','keepEns','nPrint','minLeaf','oobVarImp'},2);

% General options for Matlab ensemble construction.
% The random stream for each training iteration is derived from 'StreamSeed'
% and the training round number.
params.train.ensemble.statset = struct( 'streamSeed',  123 );



%% 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %








