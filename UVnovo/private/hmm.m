function [alpha, beta_mc, paramsOut] = hmm(probSpec, transMat, pmass_n, tm_weight)
% HMM for UVnovo peptide sequencing.
% @TODO documentation. Clean up this func. Better var names.
% probSpec: 'probability' at each integer position along a spectrum
% transMat: struct defines HMM transitions, frequencies, & n-/c-term conditions.
% pmass_n: sum of aa nominal masses composing peptide (without h2o, proton)
% paramsIn: additional parameters. (For now it's just transMat weighting).

% Parameters
if ~exist('tm_weight','var') || isempty (tm_weight)
	% transition matrix weight
	tm_weight = 0.5;
end


% Valid HMM transitions.
% These are the unique nominal residue masses that may appear in a sequence.
aamass_n = transMat.masses;
% N-/C-terminal nominal mass.
ncderiv_n = transMat.ncderiv;
% Number of possible HMM transitions.
n_aas = numel(aamass_n);

% Consts to help with preventing underflow.
BIG = 1e20;
BIGI = 1/BIG;

% This serves to normalize aa n/c-term and pair probabilities.
tm_weightConst = (1-tm_weight)/n_aas;


ssj = zeros(pmass_n,1);
ssj(probSpec(:,1)) = probSpec(:,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute forward markov chain

f1_masses = ncderiv_n(1) + aamass_n;
adjusted_f1 = ssj(f1_masses);

% Weighted n-term aas.
tm = (transMat.start*tm_weight + tm_weightConst);

statej = [f1_masses, tm.*adjusted_f1, eye(n_aas)];

maxsteps = floor(pmass_n/min(aamass_n));

alpha(1, 1:maxsteps) = struct('s', zeros(0, n_aas + 2) );
alpha(1).s = statej;

normCounter = zeros(maxsteps, 1);
j = 1;
while sum(statej(:,2)) < BIGI
	statej(:,2) = statej(:,2) * BIG;
	normCounter(j) = normCounter(j)+1;
end

mass_onpath = state_masses(transMat, pmass_n);

t = find(mass_onpath,1,'last');
massHighLims = [t, t-aamass_n(1)];
for j = 2:maxsteps
	statei = alpha(j-1).s;
	n_prev = size(statei,1);
	n_paths = n_aas;
	statej = zeros(n_prev*n_paths, n_aas+2);
	
	for n = 1:nnz(statei(:,1) <= massHighLims(2))
		
		mass_i = statei(n,1);   % mass of state i
		p_i = statei(n,2);      % probability of state i
		aa_i = statei(n,3:end); % aa contributions leading to state i
		
		mass_j = mass_i + aamass_n;
        if mass_j(end) > massHighLims(1)  % remove too-large AA masses
            mass_j(mass_j > massHighLims(1)) = [];
            n_paths = size(mass_j,1);
        end
		peaks = ssj(mass_j);
		
		t_tm = (aa_i * transMat.forward(1:n_paths, :)' )';
		tm = t_tm*tm_weight + tm_weightConst; % weight transition matrix
		
		p_j = p_i * tm .* peaks;    % prob of states j
		
		statej(n_aas*(n-1)+1:(n_aas*(n-1)+n_paths),:) = ...
			[mass_j, p_j, eye(n_paths), zeros(n_paths,(n_aas-n_paths))];
	end
	
	statej(statej(:,1)==0,:) = [];
	% Remove forward states that can't be reached backward.
	statej(~mass_onpath(statej(:,1)), :) = [];
	
	n_statej = size(statej,1);
	if n_statej == 0
		nMaxAlphaStates = j-1;
		break
	end
	while sum(statej(:,2)) < BIGI
		statej(:,2) = statej(:,2) * BIG;
		normCounter(j) = normCounter(j)+1;
	end
	
	[statej_mass, ~, ic] = unique(statej(:,1));
	n_jmass = size(statej_mass,1);
	
	% Combine state probabilities at time j with equal masses.
	tot_probs = accumarray(ic, statej(:,2));
	
	alpha(j).s = [statej_mass, tot_probs, zeros(n_jmass, n_aas)];
	
	% Normalized contribution of aa to each state.
	inds_m = accumarray(ic, 1:n_statej, [n_jmass,1], @(x){x});
	for m = 1:n_jmass
		alpha(j).s(m,3:n_aas+2) = ...
			(statej(inds_m{m},2)/tot_probs(m))' * statej(inds_m{m},3:n_aas+2);
	end
end
anormCounter = cumsum(normCounter);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute backward markov chain.

f1_masses = ncderiv_n(2) + aamass_n;
ssjB = flipud(ssj(1:pmass_n-1)); % Fliped so that it can be indexed in forward direction !
adjusted_f1B = ssjB(f1_masses);

tm = (transMat.end*tm_weight + tm_weightConst);  % Weighted c-term aas.

% Nomenclature of state i and j will be same as forward: j follows i going backward.
statej = [f1_masses, tm.*adjusted_f1B, eye(n_aas)];

beta_mc(1,1:maxsteps) = struct('s',zeros(0,n_aas+2));
beta_mc(1).s = statej;

normCounter = zeros(maxsteps,1);
j = 1;
while sum(statej(:,2)) < BIGI
	if sum(statej(:,2)) == 0
		error('hmm:ZeroSumProb',...
			'sum(states) cannot be zero. Check that probSpec was constructed properly.')
	end
	statej(:,2) = statej(:,2) * BIG;
	normCounter(j) = normCounter(j)+1;
end

mass_onpathB = mass_onpath(end-1:-1:1);

t = find(mass_onpathB,1,'last');
massHighLimsB = [t, t-aamass_n(1)];

for j = 2:maxsteps
	
	statei = beta_mc(j-1).s;
	n_prev = size(statei,1);
	n_paths = n_aas;
	statej = zeros(n_prev*n_paths, n_aas+2);
	
	for n = 1:nnz(statei(:,1) <= massHighLimsB(2))
		
		mass_i = statei(n,1);   % mass of state i
		p_i = statei(n,2);      % probability of state i
		aa_i = statei(n,3:end); % aa contributions leading to state i
		mass_j = mass_i+aamass_n;
		if mass_j(end) > massHighLimsB(1)  % remove too-large AA masses
			mass_j(mass_j > massHighLimsB(1)) = [];
			n_paths = size(mass_j,1);
		end
		peaks = ssjB(mass_j);
		t_tm = (aa_i * transMat.backward(1:n_paths, :)' )';
		tm = t_tm*tm_weight + tm_weightConst;
		
		p_j = p_i * tm .* peaks;  % prob of states j
		statej(n_aas*(n-1)+1:(n_aas*(n-1)+n_paths),:) = ...
			[mass_j, p_j, eye(n_paths), zeros(n_paths,(n_aas-n_paths))];
	end
	statej(statej(:,1)==0,:) = [];
	% Remove forward states that can't be reached backward.
	statej( ~mass_onpathB(statej(:,1)), :) = [];
	
	n_statej = size(statej,1);
	if n_statej == 0
		nMaxBetaStates = j-1;
		break
	end
	while sum(statej(:,2)) < BIGI
		statej(:,2) = statej(:,2) * BIG;
		normCounter(j) = normCounter(j)+1;
	end
	
	[statej_mass, ~, ic] = unique(statej(:,1));
	n_jmass = size(statej_mass,1);
	% Combine state probabilities at time j with equal masses.
	tot_probs = accumarray(ic, statej(:,2));
	
	beta_mc(j).s = [statej_mass, tot_probs, zeros(n_jmass, n_aas)];
	
	% Normalized contribution of aa to each state.
	inds_m = accumarray(ic, 1:n_statej, [n_jmass,1], @(x){x});
	for m = 1:n_jmass
		beta_mc(j).s(m,3:n_aas+2) = ...
			(statej(inds_m{m},2)/tot_probs(m))' * statej(inds_m{m},3:n_aas+2);
	end
end

for j = 1:nMaxBetaStates
	beta_mc(j).s(:,1) = pmass_n-beta_mc(j).s(:,1);
	beta_mc(j).s = beta_mc(j).s(end:-1:1,:);
end

bnormCounter = cumsum(normCounter);

% Convert alpha, beta to log probs.
for j = 1:nMaxBetaStates
	beta_mc(j).s(:,2) = log(beta_mc(j).s(:,2)) + log(BIGI)*bnormCounter(j);
end
for j = 1:nMaxAlphaStates
	alpha(j).s(:,2) = log(alpha(j).s(:,2)) + log(BIGI)*anormCounter(j);
end

if nargout > 2
	paramsOut.tm_weight = tm_weight;
	paramsOut.pmass_n = pmass_n;
	paramsOut.maxsteps = maxsteps;
end

end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mass_onpath = state_masses(TM, pmass_n)
	% Compute logical vector of all potential mass nodes for a given pep mass.
	% @TODO No need to call for every spectrum. Compute only once for each mass.
	massTransitions = TM.masses;
	
	maxstate = floor(pmass_n/min(massTransitions));
	
	massj = massTransitions;
	masses = massj;
	for j = 1:maxstate
		[xx, yy] = meshgrid(massj, massTransitions);
		massx = unique(sum(cat(3, xx, yy),3));
		massj = massx(~ismember(massx,masses));
		massj(massj>pmass_n) = [];
		masses = [masses; massj];
	end
	
	mForward = masses + TM.ncderiv(1);
	mBackward = -masses - TM.ncderiv(2);
	
	mass_onpath = false(pmass_n,1);
	mass_onpath(mForward(ismember(mForward, pmass_n+mBackward))) = true;
end
