function dredpreds = derivedPredictors(nodes, paramsIn)
% derivedPredictors: Create secondary predictors for RF input.
%	Return cellstr of valid predictor names when no input is provided.
% 
% INPUT
%	nodes <struct [n x 1]> Peptide interior nodes and peptide mass.
% 		.fwd <float [m(n) x 1]> peptide n-term nominal mass nodes
%		.pmass_n <float [n x 1]> nominal mass of pep residues. That is:
%			round((peptide_mass - mH2O - mProt)/unit_g)
% 
%	paramsIn <struct> input parameters
%		.appendPredictors <cellstr {1 x p}> names of predictors to return
%			valid preds: 'ionMassFrac', 'centiBinNterm', 'centiBinCterm'
%			Predictors are described below in paramsDef definition.
%			Order is maintained for output.
% 
% OUTPUT
%	dredpreds <cell [n x p] of float [m(n) x 1]> Predictor variables

paramsDef = struct( ...
	'appendPredictors', {{
		'ionMassFrac';     % ion mass / precursor mass
		'centiBinNterm';   % 100 Da window containing the N-term mass
		'centiBinCterm';   % 100 Da window containing the C-term mass
		}} ...
	);
if ~exist('paramsIn','var'), paramsIn = []; end
params = initParams(paramsDef, paramsIn);

%%%

if ~exist('nodes','var')
	dredpreds = params.appendPredictors;
	return
end

%%%
npeps = numel(nodes);

nAppendPreds = numel(params.appendPredictors);

dredpreds = cell(npeps, nAppendPreds);
for i = 1:npeps
	nodeMasses = nodes(i).fwd;
	pmass_n = nodes(i).pmass_n;
	for n = 1:nAppendPreds;
		predName = params.appendPredictors{n};
		% This has less memory overhead than defining a anonymouse function
		% for each predictor.
		switch predName
			case 'ionMassFrac'
				dredpreds{i,n} = nodeMasses ./ pmass_n;
			case 'centiBinNterm'
				dredpreds{i,n} = ceil(0.01*nodeMasses);
			case 'centiBinCterm'
				dredpreds{i,n} = ceil(0.01*(pmass_n - nodeMasses));
			otherwise
				warning(exceptionID('unknownPredVarName'), ...
					'Unrecognized predictor var name %s', predName)
				dredpreds{i,n} = nan;
		end
	end
end



%% OLD more flexible framework outline 
% % % ARGS
% % %	preds: <cell [m x 2]> {[pred name <str>,<funcHandle>];...}
% % %		Functions must be defined in terms of query mass:
% % %			'ionMasses'
% % %		and the primary spectra variables:
% % %			'pmass', 'pcharge'
% % %		and/or vars in predVecs: ?????????? useful but hard to implement cleanly
% % %		and/or predVars preceding current one that's being evaluated.
% % % % %		and/or the sample variables:  % these should be part of anon func
% % % % %			'AAs' (as from MSaalist), 
% % %	specVars: <struct>
% % %		.pmass   <array [n x 1]> precursor mass
% % %		.pcharge <array [n x 1]> precursor charge
% % %		
