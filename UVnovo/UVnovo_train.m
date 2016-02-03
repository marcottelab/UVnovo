% WORK IN PROGRESS. This is mostly an outline of a code I'm porting over.

function varargout = UVnovo_train(Meta)
% UVNOVO_TRAIN learns peptide fragmentation characteristics from known spectra
%	and constructs random forests for later interpretation of unknown spectra.
%	ms2 spectra
% 
% @TODO add documentation.
% 
% See also UVNOVO, UVNOVO_DENOVO.


% Load training data (spectra with corresponding PSM sequences).
[msData, psmData, mapi] = import_data(Meta);

% Get nominal fragment masses (nodes) for theoretical fragment ions.
nomNodes = getNominalNodes({psmData.scans.seqAnno});

% Prepare feature vectors for training observations

% Feature vectors for true observations.
fvPos = createFeatureVecs();

% Negative-example training vectors.

% Iteratively train Random Forests.
ensParams = Meta.params.train.ensemble.iters; % Ensemble params for each round.
nIterations = numel(ensParams);

rfEnsembles = struct('ens',[], 'shifts',cell(nIterations,1));
featureImportances = cell(nIterations,1);

for nIter = 1:nIterations
	% train
end

% Save RF.

varargout = {};
end %MAIN

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function featureVecs = createFeatureVecs()
	% Create vector of fragment ion feature vars for each position on each scan.
end

function dredpreds = derivedPredictors()
	% Create secondary predictors for RF input.
	% 	(see H:\Programming\matlab\REPOS\uvnovo\dev\derivedPredictors.m)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Complete functions.

function nomNodes = getNominalNodes(seqs)
	% Get nominal fragment masses (nodes) for theoretical fragment ions.
	% SEQS <cellstr> annotated peptide sequences.
	
	% Struct of amino acid, n/c-term, and elemental nominal masses.
	aas_nom = MSaalist('masstype','nominal');
	
	[uniseqs, ~, uni2allseqs] = unique(seqs);
	nUni = numel(uniseqs);
	
	nodesFwd = cell(nUni, 1);
	nodesRev = cell(nUni, 1);
	
	% Calculate theoretical ion series masses from annotated peptide sequences.
	for i = 1:nUni
		ss = MSsynthspec(uniseqs{i}, aas_nom, 'no');
		nodesFwd{i} =        ss(1:end-1, 1);
		nodesRev{i} = flipud(ss(1:end-1, 2));
	end
	
	[nomNodes(1:nUni).fwd] = nodesFwd{uni2allseqs};
	[nomNodes(1:nUni).rev] = nodesFwd{uni2allseqs};
end

function [msData, psmData, mapi] = import_data(Meta)
	% Load spectra and corresponding PSM sequences.
	% This can either be from a MS2 and psm file pair or a training .mat file.
	
	paths = Meta.paths;
	
	filetypes = fieldnames(paths);
	
	if all(ismember({'spectra', 'psms'}, filetypes))
		% Import from MS2 and psm files.
		msData = import_spectra(paths.spectra.path, Meta.params.pre.import_spectra);
		psmData = import_psms(paths.psms.path, 'verbosity', Meta.params.UVnovo.verbosity);
		
	elseif ismember('train', filetypes)
		% Load spectra and psm data from file created by UVnovo_partition.
		s = io.loadSer(paths.train.path, 'asStruct');
		msData = s.msData;
		psmData = s.psmData;
		%partitionMeta = s.Meta; % @TODO check for discrepancies with Meta that could cause a problem.
		clear s
	end
	
	% Index between psmData.scans & msData.scans.
	mapi = cross_index(psmData, msData);
	
end












