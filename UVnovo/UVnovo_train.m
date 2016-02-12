function varargout = UVnovo_train(Meta, varargin)
% UVNOVO_TRAIN learns peptide fragmentation characteristics from known spectra
%	and constructs random forests for later interpretation of unknown spectra.
%	ms2 spectra
% 
% @TODO add documentation.
% 
% See also UVNOVO, UVNOVO_DENOVO.


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% For development & debugging purposes. ARG1 may be the name of a subfunction
% instead of Meta. In this case, only that subfunction is called and any outputs
% are immediately returned.
if ischar(Meta) && strcmp(which(Meta), which(mfilename))
	fh = str2func(Meta);
	[varargout{1:nargout}] = fh(varargin{:});
	return
end


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Load training data (spectra with corresponding PSM sequences).
[msData, psmData, trainingDataFiles] = import_data(Meta);
fprintf(1, 'File loaded:%s', sprintf('\t%s\n', trainingDataFiles{:}))
fprintf(1, 'Preparing data for random forest construction ... ')

% Get nominal fragment masses (nodes) for theoretical fragment ions.
[uniseqs, ~, seqs2uni] = unique({psmData.scans.seqAnno}');
nodes = getSeqNodes(uniseqs);

% Index between elements of psmData.scans, msData.scans, and nodes.
% MAPI.(A).(B): index in B for each element in A.
mapi = struct;
% @TODO extract to cross_index(). % mapi = cross_index(psmData, msData);
psm_nscans = [psmData.scans.nscan]';
ms_nscans = [msData.scans.nscan]';
[~, mapi.psmData.msData] = ismember(psm_nscans, ms_nscans);
[~, mapi.msData.psmData] = ismember(ms_nscans, psm_nscans);
% Indexing to nodes.
mapi.psmData.nodes = seqs2uni;
mapi.msData.nodes = mapi.psmData.nodes(mapi.msData.psmData);


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Prepare feature predictor vectors for Random Forest training.

% Get initial set of predictor names.
predNames = initPredictorNames(Meta);

scans = msData.scans;
scan2node = mapi.msData.nodes;
pvparams = Meta.params.predVecs;

% Predictor vectors-- True observations.
pvPosData = createPredVecs(scans, nodes, scan2node, predNames, pvparams);

% Randomly shift nodes.
% @TODO The shifts are inconsistent with previous program!
nodesNeg = rshiftNodes(Meta, nodes);
% Predictor vectors-- Negative training examples.
pvNegData = createPredVecs(scans, nodesNeg, scan2node, predNames, pvparams);


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Make data (X, Y, W) for training TreeBagger() ensemble.

% Index in psmData for each predVec observation.
pvPos2psmData = mapi.msData.psmData( pvPosData.mapi.predVecs_expand.scans );
pvNeg2psmData = mapi.msData.psmData( pvNegData.mapi.predVecs_expand.scans );

% Downweight observations from more abundant peptides during ensemble training.
% Use sqrt(counts) rather than total counts because higher count peps likely
% generate nicer spectra & are on average higher confidence.
psmCounts = [psmData.scans.psmCountInSample]';
psmWeights = 1./sqrt(psmCounts);
wTrain = psmWeights([pvPos2psmData; pvNeg2psmData]);

% This is the training data. -- predictors and response vars.
xTrainFull = cat(1, pvPosData.pvecs.preds, pvNegData.pvecs.preds);
yTrain = [true(size(pvPos2psmData, 1), 1); false(size(pvNeg2psmData, 1), 1)];

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Iteratively train Random Forests and select best predictors from each round.

% Parameters for each training round.
ensParams = Meta.params.train.ensemble.iters;
nIterations = numel(ensParams);
featureImportances = cell(nIterations,1);

fprintf(1, 'Done.\nTraining random forests over %g iterations ...\n', nIterations)
for nIter = 1:nIterations
	params = ensParams(nIter);
	params.useParallel = Meta.params.UVnovo.useParallel;
	params.nPrint = Meta.params.train.ensemble.nPrint;
	% @TODO This rstream seed is inconsistent with previous program.
	% 	Old version: seed = rstreamSeed + 256*nPartition + nIter;
	params.rstreamSeed = Meta.params.train.ensemble.rstreamSeed + nIter;
	
	if params.nvars == 0 || params.nvars >= numel(predNames)
		% Use all predictors for current round of training.
		nvars = numel(predNames);
		varInds = true(1, nvars);
		params.nvars = nvars;
	else
		% Get subset of most important predictors for training.
		assert(nIter > 1, exceptionID('VarSelectionNotPossible'), ...
			'Must train using all vars before var importance is known.')
		nvars = params.nvars;
		
		[~, ia] = sort([featureImportances{nIter-1}{:,1}]);
		rankedPredictors = featureImportances{nIter-1}(ia,2);
		
		[a, ia] = ismember(rankedPredictors(1:nvars), predNames);
		assert(all(a));
		varInds = false(1, nvars);
		varInds(ia) = true;
	end
	
	xTrain = xTrainFull(:, varInds);
	params.predNames = predNames(varInds);
	
	% Train a Random Forest and construct ensemble 'ensTB' of class TreeBagger.
	fprintf(1,'Round %g. Building %d trees from %d predictors ...\n', ...
		nIter, params.ntrees, nvars); tic
	ensTB = trainRF(xTrain, yTrain, wTrain, params);
	fprintf(1,'\tRF constructed in %f seconds.\n', toc)
	
	% Calculate and plot out-of-bag error vs. number of trees grown.
	% The oob error calculation (part of Matlab TreeBagger class) is slow.
	if params.plotOOBError
		oobErrorVals = plotOOBError(ensTB);
		fprintf(1,'\tOOB error: %f\n', oobErrorVals(end))
	end
	
	% Get predictor importance metrics and ranks.
	if strcmpi(params.oobVarImp, 'on')
		scoreMetrics = {'OOBPermutedVarDeltaError', 'DeltaCritDecisionSplit'};
	else
		scoreMetrics = {'DeltaCritDecisionSplit'};
	end
	featureImportances{nIter} = getPredImportance(ensTB, scoreMetrics);	
end

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Save the final Random Forest.

[outdir, t_fn] = fileparts(trainingDataFiles{1});
fn = [t_fn, '_RF_', datestr(now, 'yyyymmdd') '.mat'];
fnEns = io.genUniquePath( fullfile( outdir, fn));

% % Compact the final ensemble.
% Ens = compact(ensTB);
Ens = ensTB;

io.saveSer(fnEns, 'Ens', 'featureImportances')
% @TODO Save the predictor importance tables in text format.

fprintf(1, 'UVnovo random forest training complete. Model saved to:\n\t%s', fnEns)

varargout = {};

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Placeholder --- Examples for calling useful UVnovo plotting functions.
function plotting_placeholder()
	% Plots for observing UVnovo_train data and results.
	
	% Plot RF predictors. Good for visualizing features across observations.
	plotPredictorVecs(predVec_mat, predNames, 'subsample', 0.1, 'clim', [-1,1])
	
	% Plot RF predictor importance. Custom tooltips help with inspection.
	plotPredictorImportance(ensTB, 'score', 'DeltaCritDecisionSplit')
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Subfunctions


function [msData, psmData, files] = import_data(Meta)
	% Load spectra and corresponding PSM sequences.
	% This can either be from a MS2 and psm file pair or a training .mat file.
	
	paths = Meta.paths;
	
	filetypes = fieldnames(paths);
	
	if ismember('train', filetypes)
		% Load spectra and psm data from file created by UVnovo_partition.
		files = {paths.train.path};
		s = io.loadSer(files{1}, 'asStruct');
		msData = s.msData;
		psmData = s.psmData;
		% @TODO check for discrepancies with Meta that could cause a problem.
		%partitionMeta = s.Meta;
		clear s
		
	elseif all(ismember({'spectra', 'psms'}, filetypes))
		% Import from MS2 and psm files.
		files = {paths.spectra.path; paths.psms.path};
		msData = import_spectra(files{1}, Meta.params.pre.import_spectra);
		psmData = import_psms(files{2}, 'verbosity', Meta.params.UVnovo.verbosity);
	end
	
end


function nodes = getSeqNodes(seqs)
	% Get nominal fragment masses (nodes) for peptide sequences.
	% Nodes (forward and reverse) represent the mass ladders of the peptide
	% amino acid residue series in the n- and c-terminal directions. These are
	% nominal masses (count of neutrons + protons).
	% 
	% SEQS <cellstr {n x 1}> annotated peptide sequences.
	% 
	% NODES <struct [n x 1]> peptide fragment mass nodes.
	%	.pmass_n <int> mass of all peptide residues (nominal pepmass - 19 Da)
	%	.fwd <array> forward nodes
	%	.rev <array> reverse nodes (pmass_n - nodes.fwd)

	% Struct of amino acid, n/c-term, and elemental nominal masses.
	aas_nom = MSaalist('masstype','nominal');
	
	nseqs = numel(seqs);
	
	pmasses = zeros(nseqs, 1);
	nodesFwd = cell(nseqs, 1);
	nodesRev = cell(nseqs, 1);
	
	% Calculate theoretical ion series masses from annotated peptide sequences.
	for i = 1:nseqs
		ss = MSsynthspec(seqs{i}, aas_nom, 'no');
		pmasses(i) = ss(end, 1);
		nodesFwd{i} = ss(1:end-1, 1);
		nodesRev{i} = flipud(ss(1:end-1, 2));
	end
	
	nodes = struct( ...
		'pmass_n', num2cell(pmasses), ...
		'fwd', nodesFwd, ...
		'rev', nodesRev  ...
		);
end


function predNames = initPredictorNames(Meta)
	% Create array of predictor names. The names specify the set of predictors
	% that are calculated and returned in createPredVecs().
	w = Meta.params.predVecs.featuresMassWin;
	offsets = w(1):w(2);
	mzOffsets_init = struct( ...
		'nterm', {offsets, offsets/2}, ... charge 1 & charge 2 feature offsets [Da]
		'cterm', {offsets, offsets/2}  );

	[~,offsetPredNames] = convertPredNames_mzOffsets(mzOffsets_init);
	derivedPredNames = derivedPredictors();

	predNames = [offsetPredNames(:); derivedPredNames(:)]';
end


function nodesShifted = rshiftNodes(Meta, nodes)
	% Shift true fragment node masses each by a random integer value.
	
	rs = RandStream('twister','Seed', Meta.params.predVecs.falseObs.rstreamSeed);
	
	% Make distribution from which to randomly draw shift for each ion.
	% Possible shifts from which to draw:
	w = Meta.params.predVecs.featuresMassWin;
	shiftPool = w(1):w(2);
	% Disallow +0 m/z shift.
	ignoreShifts = [0];
	shiftPool(ismember(shiftPool,ignoreShifts)) = [];
	% Weights (currently uniform) define proclivity to for certain mass shift.
	% From basic testing, non-uniform weighting schemes did not help final RFs.
	weights = ones(size(shiftPool));
	
	% Partition interval (0,1) into bins sized relative to shift weights & draw
	% uniformly to sample from the weight distribution.
	shiftCDF = cumsum(weights(:)./sum(weights))';
	sspec_sizes = cellfun('size', {nodes.fwd}, 1);
	n = sum(sspec_sizes);
	x = [rand(rs,1,n), shiftCDF];
	[~,ia] = sort(x);
	t = ia;
	t(ia) = cumsum( accumarray(find(ia>n)', 1, size(x')) ) + 1;
	shifts = shiftPool(t(1:n))';
	
	nodesShifted = nodes; % Preallocate. Values will be replaced.
	ilo = 1;
	for i = 1:numel(nodes)
		ihi = ilo+sspec_sizes(i)-1;
		nodesShifted(i).shifts = shifts(ilo:ihi);
		nodesShifted(i).fwd = nodes(i).fwd + shifts(ilo:ihi);
		nodesShifted(i).rev = nodes(i).rev - shifts(ilo:ihi);
		ilo = ihi+1;
	end
end


function ensTB = trainRF(X, Y, W, params)
	% Train Random Forest ensemble.
	
	statsetOpts = statset( ...
		'useParallel', params.useParallel, ...
		'streams', RandStream('multFibonacci','Seed', params.rstreamSeed), ...
		'useSubstreams', params.useParallel ...
		);
	
	ensVarargs = {
		'names',     params.predNames;
		'weights',   W;
		'options',   statsetOpts;
		'oobVarImp', params.oobVarImp;
		'nPrint',    params.nPrint;
		'minLeaf',   params.minLeaf;
		'method',   'classification';
		'oobPred',  'on'; % needed for oobError(ensTB). No apparent speed penalty.
		% 'cat',    []; % categorical var indices (none at present)
		}';
	
	ensTB = TreeBagger(params.ntrees, X, Y, ensVarargs{:});
end


function predImportanceTable = getPredImportance(ensTB, scoreMetrics)
	% Get importance metric for each var in TreeBagger ensemble 'ensTB'.
	% ENSTB <TreeBagger or CompactTreeBagger>
	% SCOREMETRICS <cellstr> built in metric(s) of predictor importance:
	%	{ 'OOBPermutedVarDeltaError', 'DeltaCritDecisionSplit',  
	%	  'OOBPermutedVarDeltaMeanMargin', 'OOBPermutedVarCountRaiseMargin' }
	%	The last two metrics are nearly identical to first, from what I've seen.
	%	See doc stats.TreeBagger for description of each.
	% 
	% OUTPUT predImportanceTable <cell array {nVars x (2 + 2*nMetrics)>
	%		Cell of predictor rank, name, and column pairs for each score
	%		metric. Each pair is metric rank & metric score.
	
	assert( isa(ensTB, 'TreeBagger') || isa(ensTB, 'CompactTreeBagger') )
	
	varNames = ensTB.VarNames;
	nvars = numel(varNames);
	nmetrics = numel(scoreMetrics);
	
	scores = zeros( nvars, 2*nmetrics);
	for i = 1:nmetrics
		% Predictor scores.
		scores(:, 2*i) = ensTB.(scoreMetrics(i))';
		% Ranks predictors.
		[~, ia] = sort( scores(:, 2*i), 'descend' );
		scores(ia, 2*i-1) = 1:nvars;
	end
	
	% Overall predictor rank is geometric mean of metric-wise ranks.
	ranks = power( prod( scores(:, 1:2:end), 2), 1/nmetrics);
	
	predImportanceTable = [num2cell(ranks), varNames', num2cell(scores)];
end


