function varargout = UVnovo_denovo(Meta, varargin)
% UVNOVO_DENOVO
% 
% @TODO add documentation.
% 
% See also UVNOVO, UVNOVO_TRAIN.


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
% Load spectra and UVnovo trained models.
fprintf(1, 'Importing spectra and UVnovo models ... ')
% AAs_HMM: struct of amino acid, n/c-term, and elemental nominal masses for HMM.
[msData, Ens, TransMat, AAs_HMM] = import_data(Meta);
fprintf(1, 'Done.\n')


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Prepare feature vectors for Random Forest predictions.

scans = msData.scans;
% Nominal mass of residues composing the peptide.
pmass_n = round(([scans.pmass]' - CONSTS.mProton - CONSTS.mH2O)/CONSTS.unit_g);

% Get valid mass nodes.
% These are the possible fragmentation sites for each spectrum. Applying the
% random forest at each of these positions yields the fragment site predictions.
[pmass_n_uni, ~, scan2node] = unique(pmass_n);
nodes = getMassNodes(pmass_n_uni, AAs_HMM);

% Make the predictor vectors (features) for each node.
varNames = Ens.VarNames;
pvparams = Meta.params.predVecs;
predVecs = createPredVecs(scans, nodes, scan2node, varNames, pvparams);


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Apply RF for peptide spectrum fragmentation site prediction.

% Select scans for RF and HMM processing. We use all for now. This will make
% it easier for processing singles or specific subsets.
indPredictSpec = 1:numel(scans);
pv = {predVecs.pvecs(indPredictSpec).preds}';

% Get predictions.
rfscores = predictFragSites(Ens, pv);



% HMM de novo sequencing...
% 
% 
% 



varargout = {};

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Subfunctions


function rfscores = predictFragSites(Ens, pv)
	% Apply RF for fragmentation site predictions.
	% Ens: <ensemble object> Random Forest used for predictions.
	% pv: <cell array {n x 1}> each cell contains a matrix of predictor vectors
	%	for a single spectrum [nObs x nVars].
	
	% RF predictions are MUCH faster when called fewer times on larger data.
	% Aggregate predictor vectors and split into large chunks for processing.
	pVecs_expand = cat(1, pv{:});
	npvecs = size(pVecs_expand,1);
	
	% Maximum number of points to predict at a time:
	maxChunkSize = 100000;
	nchunks = ceil(npvecs/maxChunkSize);
	
	% Get [low, high] indices for each predict chunk:
	t = linspace(0, npvecs, nchunks + 1)';
	chunks = round([t(1:end-1) + 1, t(2:end)]);
	
	rfscores_expand = zeros(npvecs, 1);
	fprintf(1, 'Random forest fragmentation site prediction ')
	tic
	for i = 1:nchunks
		ilo = chunks(i,1);
		ihi = chunks(i,2);
		[~, s] = predict(Ens, pVecs_expand(ilo:ihi, :));
		rfscores_expand(ilo:ihi) = s(:,2);
		
		fprintf(1, '.')
	end
	t = toc;
	fprintf(1, ['\nEnsemble predictions complete in %f seconds,\n'...
		'\t%f points/second\n'], t, ihi/t);
	
	% Index back to individual pv cells from concatenated predictor vectors.
	npeaks = cellfun('size', pv, 1);
	pvec2scan = cumsum(accumarray( cumsum([1; npeaks(1:end-1)]), 1, [npvecs,1]));
	
	% Allocate predictions back into cells of correct size.
	rfscores = accumarray(pvec2scan, rfscores_expand, [], @(x){x});
end


function nodes = getMassNodes(pmass_n, AAs_HMM)
	% Get all potential fragment masses (nodes) for peptides of given mass.
	% Nodes (forward and reverse) represent the nominal masses in the n- and
	% c-terminal directions.
	% 
	% PMASS_N <array [n x 1]> nominal mass of peptide aas (nom pepmass - 19 Da).
	% 
	% NODES <struct [n x 1]> possible peptide fragment mass nodes.
	%	.pmass_n <int>
	%	.fwd <array> forward nodes
	%	.rev <array> reverse nodes (pmass_n - nodes.fwd)
	
	nmasses = numel(pmass_n);
	
	minAAintmass = min(AAs_HMM.aamass);
	minNode = minAAintmass + AAs_HMM.ncderiv(1);
	maxNodes = pmass_n - minAAintmass - AAs_HMM.ncderiv(2);
	
	nodesFwd = cell(nmasses,1);
	nodesRev = cell(nmasses,1);
	for n = 1:nmasses
		nodesFwd{n} = (minNode:maxNodes(n))';
		nodesRev{n} = pmass_n(n) - nodesFwd{n};
	end
	
	nodes = struct( ...
		'pmass_n', num2cell(pmass_n), ...
		'fwd', nodesFwd, ...
		'rev', nodesRev  ...
		);
	
end


function [msData, Ens, TransMat, AAs_HMM, filesIn] = import_data(Meta)
	% Load spectra, random forest, and initialize aa mass transition matrix.
	
	paths = Meta.paths;
	filetypes = fieldnames(paths);
	filesIn = cell(3,1);
	if ismember('test', filetypes)
		% Load spectra from test file created by UVnovo_partition, or any
		% serialized MAT file that alreday contains a 'msData' var.
		fn = paths.test.path;
		filesIn{1} = fn;
		s = io.loadSer(fn, 'asStruct');
		if ~isfield(s, 'msData')
			error('UVnovo_denovo:import_data:msDataNotInMatFile', ...
				'File must contain the msData spectra variable.\n\t%s', fn)
		end
		msData = s.msData;
		
	elseif ismember('spectra', filetypes)
		% Import from MS2 file.
		filesIn{1} = paths.spectra.path;
		msData = import_spectra(filesIn{1}, Meta.params.pre.import_spectra);
	end
	
	% Import random forest.
	filesIn{2} = paths.rf.path;
	s = io.loadSer(filesIn{2}, 'asStruct');
	Ens = s.Ens;
	clear s
	
	% Init residue mass transition matrix.
	filesIn{3} = paths.aamodel.path;
	[TransMat, AAs_HMM] = massTransMat(filesIn{3}, Meta.params.AAs);
	
end

