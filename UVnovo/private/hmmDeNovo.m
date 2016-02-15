function denovoData = hmmDeNovo(rfpeaks, TM, paramsIn)
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
	'peakScoreRecordThresh', 1e-4, ... Remove HMM peaks below x.
	'minQuickScore_to_proceed', -5 ... Skip spectra with low quick scores.
	);
if ~exist('paramsIn','var'), paramsIn = []; end
params = initParams(paramsDef, paramsIn);


% @TODO get rid of the following fp (fragment prediction) templates.
t ={
	[],                'predSeq';
	[],                'peakScores';
	[],                'nodeMasses';
	[],                'nodeScores';
	0,                 'avgNodeScore';
	0,                 'minNodeScore';
% 	[],                'r_nodeScores';
% 	0,                 'r_avgNodeScore';
% 	0,                 'r_minNodeScore';
	};
fpTemplate = cell2struct( t(:,1), t(:,2), 1);
t ={
	[],                'peakScores';
	[],                'nodeScores';
	0,                 'avgNodeScore';
	0,                 'avgNodeScore_normalized';
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

for m = 1:numel(rfpeaks)
	pmass_n = rfpeaks(m).pmass_n;
	probSpec = [rfpeaks(m).mass, (rfpeaks(m).score - 0.5)*params.probRegularize + .5];
	specForHMM = [probSpec(:,1), probSpec(:,2)./(1-probSpec(:,2))];
	
	[alpha, beta_mc] = hmm(specForHMM, TM, pmass_n, params.tm_weight);
	
	
	% Compile stats on predictions at different numbers of fragments.
	% Get bounds for realizeable sequences.
	% @TODO calculate this once for each pep mass.
	% Lower length bound:
	nfInit = ceil( (pmass_n - sum(TM.ncderiv) - 1) / max(TM.masses));
	while ~any( ismember(alpha(1).s(:,1), beta_mc(nfInit-1).s(:,1)) )
		nfInit = nfInit+1;
	end
	% Higher length bound:
	nfFinal = numel(beta_mc);
	while ~any(ismember(alpha(1).s(:,1), beta_mc(nfFinal-1).s(:,1)))
		nfFinal = nfFinal-1;
	end
	
	% %% Quick-score.
	
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
		t = total;
		t( isnan(t) | (t < params.peakScoreRecordThresh) ) = 0;
		fpQuick(n_frags).peakScores  = sparse(t);
		fpQuick(n_frags).nodeScores  = total(state_peaks);
		fpQuick(n_frags).avgNodeScore  = mean(total(state_peaks));
	end
	t = [fpQuick(nfInit:nfFinal).avgNodeScore]';
	t = num2cell( (t - min(t)) / range(t) );
	[fpQuick(nfInit:nfFinal).avgNodeScore_normalized] = t{:};
	
	t = [fpQuick.avgNodeScore];
	quickScore = max( t(logical(t)) );
	
	denovoData(m).quickscore = fpQuick;
	
	if quickScore < params.minQuickScore_to_proceed
		fprintf(1,'spectrum quickscore: %f.  Score below threshold for sequencing.\n%d spectra completed.\n',...
			quickScore, m);
		continue
	else
		msgBuffer = sprintf('spectrum quickscore: %f\n%d spectra completed.\n',...
			quickScore, m);
	end
	
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
	for n = 1:ntops % @TODO Repetitive code could be optimized. Profile first.
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
			%	alpha_i where overlaps beta
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

		% Retaining peak scores is needed to rescoring. @TODO move the
		%	rescoring function here & make optional to retain peak scores.
		t = total_eachstate;
		t(t < params.peakScoreRecordThresh) = 0;
		fp(n_frags).peakScores = sparse(t);

		%%% Rescoring
		% @TODO move rescoring function here. Replace inefficient one below.
		if m == 1 && n == 1
			warning('Fill in the sequence rescoring stuff!!!!!')
		end
		%	Nodes that match any non-current predFrags site are removed for score recalculations.
		% % % % 			r_total_eachstate = zeros(pmass_n,n_frags-1);
		% % % % 			r_predNodeScores = zeros(n_frags-1,1);
		% % % % 			for i = 1:n_frags-1
		% % % % 				state_a = alpha(i).s(indsAlphaInBeta{n,i},:);
		% % % % 				state_b = beta_mc(n_frags-i).s(indsBetaInAlpha{n,i},:);
		% % % % 				oknodes = ~ismember(state_a(:,1), predFrags([1:i-1, i+1:end]));
		% % % %
		% % % % 				lo_sum = state_a(oknodes,2)+state_b(oknodes,2);
		% % % % 				t = exp(lo_sum-max(lo_sum));
		% % % % 				stateScores = t/sum(t);
		% % % % 				r_total_eachstate(state_a(oknodes,1),i) = stateScores;
		% % % % 				r_predNodeScores(i) = r_total_eachstate(predFrags(i),i);
		% % % % 			end
		% % % % 			fp(n_frags).r_nodeScores = r_predNodeScores;
		% % % % 			fp(n_frags).r_avgNodeScore = mean(r_predNodeScores);
		% % % % 			fp(n_frags).r_minNodeScore = min(r_predNodeScores);
	end
	
	% Get top-scoring sequecne prediction.
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
	
	fprintf(1,'%s', msgBuffer)
end

