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

% Ensemble training schedule. Iteratively train & select top predictor vars for
% each subsequent model.
% Parameters for each training round:
%	nvars: Select this many predictor vars. First step must be '0' (all vars).
%	ntrees: Number of decision trees in the ensemble.
%	minLeaf: Minimum number of observations per tree leaf.
%	oobVarImp: <'on'|'off'> Out-of-bag estimate of var importance. SLOW!
%	plotOOBError: Calculate & plot out-of-bag error vs. ntrees. It can be slow.
% User may add or remove rows from following table as desired.
params.train.ensemble.iters = cell2struct({
	...<nvars>  <ntrees>  <minLeaf>  <oobVarImp>  <plotOOBError>
		0,       150,      1,        'off',        true;    % round 1
		200,     150,      1,        'off',        true;    % round 2
		100,     150,      1,        'off',        true;    % round 3
		60,      250,      1,        'on',         true;    % round 4
		30,      400,      1,        'on',         true;    % final round
	}, ...
	{  'nvars', 'ntrees', 'minLeaf', 'oobVarImp', 'plotOOBError'}, 2);

% Print progress after for every n trees trained.
params.train.ensemble.nPrint = 4;

% Random stream seed for Matlab's ensemble training procedure.
params.train.ensemble.rstreamSeed = 123;

%% Random Forest predictor (feature) variables.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %

% Parameters controlling how ms2 peaks are normalized and converted to features:
% ms2 fragment ion mass accuracy [m/z]. Peaks within this range are binned.
params.predVecs.mzThreshold = [-.5, .5];
% Window [m/z] for local peak intensity normalization.
params.predVecs.localPeakWin = [-50, 50];

% Mass window [Da] defines which peaks, around each nominal frag site, are
% considered as features. This is only for the first training round, before the
% important predictors are identified.
params.predVecs.featuresMassWin = [-50 50];

% Random stream seed value for repeatable generation of 'random' true-negative
% training examples. If left empty, UVnovo assigns a random seed, and the
% training procedure will yield slightly different results each time it is run.
params.predVecs.falseObs.rstreamSeed = 121212;






