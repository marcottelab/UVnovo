function [denovoData, params] = hmmDeNovo(rfpeaks, TM, paramsIn)
% HMMDENOVO ...
%	...
% 
% INPUT ARGUMENTS
%	RFPEAKS <struct [n x 1]> peptide fragment site RF predictions.
%		.pmass_n <int> mass of residues composing a peptide (nom pepmass - 19Da)
% 		.mass    <array [m(n) x 1]> mass of forward nodes (peaks)
% 		.score   <array [m(n) x 1]> rf score for each peak
%	AAs <struct> Residue, n/c-term, and elemental nominal masses.
%	TM  <struct> Residue frequencies and transition matrices.
%	PARAMSIN <struct> HMM parameters, described below for 'paramsDef'.
% 
% OUTPUT
%	***: <struct [n x 1]> ***
% 
% @TODO documentation, code cleanup, refactor output var structure.


paramsDef = struct( ...
	'probRegularize', .999,... RF probability moderation [0.9 < x < 1].
	'fastFragThresh', .92, ... Quickscore cutoff for sequencing [0 <= x <= 1].
	'tm_weight', .5, ...       AA transition matrix weight [0 <= x <= 1].
	'peakScoreRecordThresh', 1e-4, ... Remove HMM peaks below x. [0 <= x <= 1]. x > 1: discard all.
	'minQuickScore_to_proceed', -5,... Skip spectra with low quick scores.
	'nodePolarizationVals', ... For each sequence length, count number of HMM nodes scored above these values.
		sort(reshape(bsxfun(@times, [1 2 5]', 10.^(-6:-1)), [], 1), 'descend') ...
	);
if ~exist('paramsIn','var'), paramsIn = []; end
params = initParams(paramsDef, paramsIn);


if params.peakScoreRecordThresh >= 0 && params.peakScoreRecordThresh <= 1
	retainPeakScores = true;
else
	retainPeakScores = false;
end

% Templates for de novo results structs.
t ={
	[],         'predSeq';
	[],         'peakScores';
	[],         'nodeMasses';
	[],         'nodeScores';
	0,          'avgNodeScore';
	0,          'minNodeScore';
	[],         'r_nodeScores';
	0,          'r_avgNodeScore';
	0,          'r_minNodeScore';
	[],         'polarization';
	};
fpTemplate = cell2struct( t(:,1), t(:,2), 1);
t ={
	[],         'nodeScores';
	0,          'avgNodeScore';
	0,          'avgNodeScore_normalized';
	};
fpQuickTemplate = cell2struct( t(:,1), t(:,2), 1);
 

% Initialize output data struct.
% @TODO Rename fields. Remove unneeded fields.
denovoData = repmat( struct( ...
		'topPreds',  struct, ...
		'fragPreds', struct, ...
		'topSeq',    '', ...
		'quickscore',struct  ...
	), size(rfpeaks) );


nodePolarizationVals = params.nodePolarizationVals;
nnpv = numel(nodePolarizationVals);

ntotal = numel(rfpeaks);

for m = 1:ntotal
	pmass_n = rfpeaks(m).pmass_n;
	probSpec = [rfpeaks(m).mass, (rfpeaks(m).score - 0.5)*params.probRegularize + .5]; %#ok<PFBNS>
	specForHMM = [probSpec(:,1), probSpec(:,2)./(1-probSpec(:,2))];
	
	[alpha, beta_mc] = hmm(specForHMM, TM, pmass_n, params.tm_weight);
	
	
	% Compile stats on predictions at different numbers of fragments.
	% Get bounds for realizeable sequences.
	% @TODO calculate this once for each pep mass.
	% Lower length bound:
	nfInit = ceil( (pmass_n - sum(TM.ncderiv) - 1) / max(TM.masses));
	while ~any( ismember(alpha(1).s(:,1), beta_mc(nfInit - 1).s(:,1) ))
		nfInit = nfInit+1;
	end
	% Higher length bound:
	nfFinal = numel(beta_mc);
	while ~any( ismember(alpha(1).s(:,1), beta_mc(nfFinal - 1).s(:,1) ))
		nfFinal = nfFinal-1;
	end
	
	
	% Quick-score.
	
	fpQuick = repmat(fpQuickTemplate, nfFinal, 1);
	% Estimate the most likely sequence lengths.
	for n_frags = nfInit:nfFinal
		state_peaks = zeros(n_frags-1, 1);
		total  = nan(pmass_n, 1);
		for i = 1:n_frags-1
			[inds_a, loc_b] = ismember(alpha(i).s(:,1), beta_mc(n_frags-i).s(:,1));
			state_a = alpha(i).s(inds_a,:);     % alpha_i where overlaps beta
			state_b = beta_mc(n_frags-i).s(loc_b(inds_a),:);
			
			lo_sum = state_a(:,2) + state_b(:,2);
			[~,imaxLO] = max(lo_sum);
			t = lo_sum;
			state_prob = t;
			
			total(state_a(:,1)) = max(total(state_a(:,1)), state_prob);
			state_peaks(i) = state_a(imaxLO,1);
		end
		% 		t = total;
		% 		t( isnan(t) | (t < params.peakScoreRecordThresh) ) = 0;
		% 		fpQuick(n_frags).peakScores  = sparse(t);
		fpQuick(n_frags).nodeScores  = total(state_peaks);
		fpQuick(n_frags).avgNodeScore  = mean(total(state_peaks));
	end
	t = [fpQuick(nfInit:nfFinal).avgNodeScore]';
	t = num2cell( (t - min(t)) / range(t) );
	[fpQuick(nfInit:nfFinal).avgNodeScore_normalized] = t{:};
	
	t = [fpQuick.avgNodeScore];
	quickScore = max( t(logical(t)) );
	
	denovoData(m).quickscore = fpQuick;
	
	msgBuffer = sprintf('%-6d Spectrum quickscore: %f', m, quickScore);
	if quickScore < params.minQuickScore_to_proceed
		fprintf(1,'%s -- Score below threshold for sequencing.\n', msgBuffer)
		continue
	end
	
	% Recalculate HMM and find best sequence (path) through HMM nodes.
	
	% Peak 'likelihoods' are transformed differently than for the quick score above.
	% @TODO Extract this magic number & figure out a good way to optimize it.
	x = .5; % x==1 is used for the quick scoring.
	specForHMM = [probSpec(:,1),  (probSpec(:,2) ./ (1 - probSpec(:,2))).^x ];
	
	[alpha, beta_mc] = hmm(specForHMM, TM, pmass_n, params.tm_weight);
	
	% Score best path for all n_frags scoring >= threshold (relative to top score).
	tops = find([fpQuick.avgNodeScore_normalized] >= params.fastFragThresh);
	ntops = numel(tops);
	indsAlphaInBeta = cell(ntops, tops(ntops)-1);
	indsBetaInAlpha = cell(ntops, tops(ntops)-1);
	for n = 1:ntops
		% @TODO Code could be optimized. Profile to see if it's a problem.
		n_frags = tops(n);
		for i = 1:n_frags-1
			[indsAlphaInBeta{n,i}, locb] = ...
				ismember(alpha(i).s(:,1), beta_mc(n_frags-i).s(:,1));
			indsBetaInAlpha{n,i} = locb(indsAlphaInBeta{n,i});
		end
	end
	
	fp = repmat(fpTemplate,	nfFinal, 1);
	for n = 1:ntops
		n_frags = tops(n);
		
		total_eachstate = zeros(pmass_n,n_frags-1);
		smass = cell(1,n_frags-1);
		for i = 1:n_frags-1
			% alpha_i where overlaps beta
			state_a = alpha(i).s(indsAlphaInBeta{n,i},:);
			state_b = beta_mc(n_frags-i).s(indsBetaInAlpha{n,i},:);
			
			smass{i} = state_a(:,1);
			
			lo_sum = state_a(:,2)+state_b(:,2);
			t = exp(lo_sum-max(lo_sum));
			stateScores = t/sum(t);
			total_eachstate(state_a(:,1),i) = stateScores;
		end
		[predFrags, predSeq, predNodeScores] = ...
			best_path(smass, total_eachstate, n_frags, pmass_n, TM);
		
		fp(n_frags).predSeq = predSeq;
		fp(n_frags).nodeScores = predNodeScores;
		fp(n_frags).nodeMasses = predFrags(1:end-1);
		fp(n_frags).avgNodeScore = mean(predNodeScores);
		fp(n_frags).minNodeScore = min(predNodeScores);
		
		% Rescoring.
		% Nodes that share a mass position with a node on the best path but
		% have a different sequence position are removed. Best-path node scores
		% are updated to reflect this.
		nodeReScores = rescorePath(predFrags(1:end-1), total_eachstate);
		fp(n_frags).r_nodeScores = nodeReScores;
		fp(n_frags).r_avgNodeScore = mean(nodeReScores);
		fp(n_frags).r_minNodeScore = min(nodeReScores);
		
		% Count of HMM nodes above certain threshold values.
		% More polarization -> node scores fall off faster. This correlates with
		% better sequence predictions.
		t = total_eachstate;
		t(predFrags,:) = 0; % ignore mass positions on best path.
		[~,ia] = sort([nodePolarizationVals; t(t > 0)], 'descend');
		ib = 1:numel(ia);
		ib(ia,1) = ib;
		fp(n_frags).polarization = ib(1:nnpv) - (1:nnpv);
		
		if retainPeakScores
			% HMM node scores. These can be used for diagnostics or
			% introspection later on.
			t = total_eachstate;
			t(t < params.peakScoreRecordThresh) = 0;
			fp(n_frags).peakScores = sparse(t);
		end
		
	end
	
	% Get top-scoring sequence prediction.
	n_frags = find([fp.avgNodeScore] == max([fp.avgNodeScore]), 1, 'last');
	predSeq = fp(n_frags).predSeq;
	
	denovoData(m).topPreds = struct( ...
		'predSeq', predSeq, ...
		'nFrags',  n_frags, ...
		'allPredSeqs', {{fp(tops).predSeq}}, ...
		'nFragsAll', tops ...
		);
	denovoData(m).fragPreds = fp;
	denovoData(m).topSeq = predSeq;
	
	fprintf(1,'%s\n', msgBuffer)
end

end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Subfunctions


function pathReScores = rescorePath(pathMasses, peakScores)
	% Nodes may have a score at more than one position along the sequence.
	% This rescoring removes the score component of out-of-position sequence
	%	nodes for each in-position sequence node.
	
	indPathNodes = sub2ind(size(peakScores), pathMasses', 1:numel(pathMasses));
	pathScores = peakScores(indPathNodes);

	% score of nodes that are on path at a different position
	t = sum(peakScores(pathMasses,:), 1) - pathScores;

	% remove score influence of the out-of-position nodes
	pathReScores = full(pathScores./(1-t))';
end

