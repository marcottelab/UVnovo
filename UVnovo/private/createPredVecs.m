function predVecs = createPredVecs(scans, nodes, scan2node, predNames, paramsIn)
% CREATEPREDVECS construct vectors of fragment ion predictor variables.
% 
% For each predefined position in a spectrum, construct the set of predictor
%	variables enumerated in predNames:
%	mzOffset predictors:
% 		Each of these features represents the tallest peak within a narrow bin
% 		offset from a specific position (node) along the peptide spectrum. The
% 		position is particular to an N- or C-terminal AA chain and ion charge.
% 		The primary position corresponds to the nominal mass of a charge 1+
% 		N-terminal chain of residues, without the extra hydrogen of N-term
% 		peptide fragments.
%		ex. 'cTerm_z1_19mz' represents y-ion charge 1+ fragments.
%	secondary predictors:
%		'ionMassFrac':   [ion_mass / precursor_mass]
%		'centiBinNterm': 100 Da window containing the N-term mass
%		'centiBinCterm': 100 Da window containing the C-term mass
% 
% INPUT
% 	scans <n x 1 struct> (only required fields are listed)
% 		.peaks
% 			.mz
% 			.rank
% 	.nodes <struct> peptide fragment mass nodes.
%		.pmass_n <double> 
% 		.fwd <array> forward nodes
% 		.rev <array> reverse nodes (pmass_n - nodes.fwd)
% 	scan2node <1 x n array> Index in nodes for each element in scans.
%	predNames <cellstr> name of predictor for each ensemble feature vector.
%	paramsIn <struct>
%		.mzThreshold  [1 x 2] ms2 m/z peak window.
%		.localPeakWin [1 x 2] size of peak intensity normalization window (Da).
% 
% OUTPUT
%	predVecs <scalar struct> predictors for ensemble training/prediction.
%		.varNames {1 x p cellstr} predictor names (should match predNames)
%		.meta.params
%				.mzThreshold
%				.localPeakWin
%		.mapi
%			.predVecs_expand.scans
%			.scans.predVecs_expand {y x 1}
%		.pvecs [n x 1]
%			.nodes [nnode(n) x 1] n-term node nominal mass (same as nodes.fwd)
%			.preds [nnode(n) x p] variable vectors for ensemble prediction
% 
% 
%	Much of this code, and that in its dependency functions, was derived from
%	other more general functions. Predictor vector construction could be much
%	more efficient, mainly by replacing mzBinVals().
% @TODO This is too slow!! Parallel processing greatly speeds mzBinVals() & is
%	used if matlabpool open.


paramsDef = struct( ... these params defined in Meta.params.predVecs
	'mzThreshold', [], ... ms2 peak mz window
	'localPeakWin', [] ... size of normalization window (Da)
	);
if ~exist('paramsIn','var'), paramsIn = []; end
params = initParams(paramsDef,paramsIn);


%% Construct mzOffset component of ensemble predictor variables (features).

% Define predictor mass offsets given predictor names.
[mzOffsets, offsetPreds, otherPreds] = convertPredNames_mzOffsets(predNames);

specToBin = [scans.peaks]';
nscans = numel(specToBin);

% mzBinVals, called below, requires that nodes be stored differently.
nodesFwdRev = cellfun(@(x,y)[x,y], {nodes.fwd}, {nodes.rev}, 'UniformOutput', false)';

mzOffsetCounts = cellfun('prodofsize',squeeze(struct2cell(mzOffsets)))';
t = num2cell(zeros(size(mzOffsetCounts)));
t(~mzOffsetCounts) = {zeros(1,0)};
zeroOffsets = cell2struct(t,{'nterm','cterm'},2);

% There will be multiple ensemble prediction vectors (observations) for each
% scan -- one for each mass position used during training or prediction.
% Create indices between the prediction vectors and scans.
t = cellfun('size', nodesFwdRev, 1);
npv = t(scan2node); % number of predictor vectors per scan
pvec2scan = cumsum( accumarray(cumsum([1;npv(1:end-1)]), 1, [sum(npv),1] ) );
scan2pvec = accumarray(pvec2scan, 1:numel(pvec2scan), [], @(x){x});


mzThresh = params.mzThreshold;
localThresh = params.localPeakWin;

% Bin peaks at each desired position and offset along all spectra.
% @TODO do this in chunks when memory starts to be an issue.
% @TODO Cell2Vec (Jan Simon file exchange) is much faster than cell2mat but
%	returns vector, not 2D matrix (easy to fix).
accumArgs = {'rank',@min,nan};
mzrankPreds = cell2mat( ...
	mzBinVals( specToBin, nodesFwdRev, scan2node, mzOffsets, mzThresh, accumArgs) ...
	);

% Get the min rank (highest intensity) peak in a sliding window along spectra.
% wide thresh -> a lot slower!! Fix this!
accumArgs = {'rank',@min,0};
localMinRanks = cell2mat( ...
	mzBinVals( specToBin, nodesFwdRev, scan2node, zeroOffsets, localThresh, accumArgs) ...
	);

% Index into localMinRanks to map back to mzranks, for peak score normalizing.
lmr2mzr = cumsum( accumarray([1,1+unique(cumsum(mzOffsetCounts(1:end-1)))]', 1, ...
	[sum(mzOffsetCounts(:)),1]))';
pmass_n = [nodes(scan2node).pmass_n];


% Normalize peak ranks by the lowest-ranked neighboring peaks. The mzOffset
% features are contained in 'mzenuf'.
%	For each peak bin:  enuf = 1 - (binRank - localMinRank + 1)/pmass_n
mzenuf = 1 - bsxfun(@rdivide, mzrankPreds - localMinRanks(:,lmr2mzr) + 1, ...
	pmass_n(pvec2scan)');
mzenuf(isnan(mzenuf)) = -1; % @TODO extract const

clear mzrankPreds localMinRanks specToBin

%% Get derived predictor vars.
derivedPreds = cell2mat( ...
	derivedPredictors(nodes(scan2node), {'appendPredictors', otherPreds})...
	);


%% Organize predictor data and return.
% Construct output struct.
predVecs = struct( ...
	'varNames', {predNames}, ...
	'meta', struct('params', params), ...
	'pvecs', struct( ...
		'nodes', {nodes(scan2node).fwd}, ...
		'preds', []) ...
	);

% Reorder predictors to match original predNames (if needed).
predNamesNew = [offsetPreds, otherPreds];
[e, ia] = ismember(predNames, predNamesNew);
assert(all(e) && length(predNames)==length(predNamesNew), ...
	exceptionID('predNameMismatch'), 'Predictors do not match original.')

% Put predictors into one big array & allocate to each spec.
preds = [mzenuf, derivedPreds];
clear mzenuf derivedPreds

if ~issorted(ia)
	for i = 1:nscans
		predVecs.pvecs(i,1).preds = preds(scan2pvec{i},ia);
	end
else
	for i = 1:nscans
		predVecs.pvecs(i,1).preds = preds(scan2pvec{i},:);
	end
end

% Map indices between scans and the concatenated prediction vectors.
mapi_pv.predVecs_expand.scans = pvec2scan;
mapi_pv.scans.predVecs_expand = scan2pvec;
predVecs.mapi = mapi_pv;

