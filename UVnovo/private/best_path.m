function [bp, phragSeq, nodeScores] = best_path(smass, sprobs, N, pmass_n, TM)
% BEST_PATH finds most likely path through nodes of HMM.
% 
%	SMASS <cell array> potential masses of each internal state.
%	SPROBS is probability at each sequence and mass position of actual peak.
%	N <int> number of residues composing sequence (# internal nodes - 1).
%	PMASS_N <int>
%	TM <struct>

aaint = TM.masses;

state = struct;
n = 1;

state(n).mass = smass{n};               % mass of each possible state n
state(n).prob = sprobs(state(n).mass,n);% probability at possible states n
state(n).leng = -log(state(n).prob);    % path length to states n
state(n).m = size(state(n).mass,1);     % number of nodes at state n
state(n).prev = ones(state(n).m,1);     % pointer to state n-1 through which shortest path to n

for n = 2:N-1
	state(n).mass = smass{n};
	state(n).prob = sprobs(state(n).mass,n);
	state(n).m = size(state(n).mass,1);
	state(n).prev = zeros(state(n).m,1);
	state(n).leng = zeros(state(n).m,1);
end

% shortstack: n, m, length, mass
shortstack = sortrows( ...
	[ones(state(1).m, 1), (1:state(1).m)', state(1).leng, state(1).mass], 3);

shortest = shortstack(1,:); % current node
n = shortest(1);
% Once n gets to N-1, shortest path is found.
while n < N-1
	shortstack(1,:) = [];
	
	% Get possible nodes(n+1,m) that aren't yet measured.
	n_hits = find( ismember( state(n+1).mass, shortest(4) + aaint ) ...
		& ~logical(state(n+1).prev) );
	state(n+1).leng(n_hits) = shortest(3) + -log(state(n+1).prob(n_hits));
	
	% Shortest path index m for state n leading to nodes in state n+1.
	state(n+1).prev(n_hits) = shortest(2);
	shortstack = sortrows([ shortstack;
		repmat(n+1, size(n_hits,1), 1), n_hits, state(n+1).leng(n_hits), ...
		state(n+1).mass(n_hits) ], 3);
	
	shortest = shortstack(1,:);
	n = shortest(1);  % current state
end

% Mass of predicted fragment ion for 0, each n, and end:
phrags = zeros(N+1,1);
phrags(1) = 0;
phrags(end) = pmass_n;

% Score for each internal node on path:
nodeScores = zeros(N-1,1);
m = shortest(2);
% Backtrack through states to get sequence.
for n = N-1:-1:1
	phrags(n+1) = state(n).mass(m);
	nodeScores(n) = state(n).prob(m);
	m = state(n).prev(m);
end

% Create string of predicted sequence.
d = diff([ TM.ncderiv(1); phrags(2:end-1); phrags(end)-TM.ncderiv(2) ]);
t_predSeq = cell(1,N);
for i = 1:numel(d)
	if TM.mass2symbol.isKey(d(i))
		v = TM.mass2symbol(d(i));
		if numel(v)>1
			t_predSeq{i} = sprintf('[%s]',v);
		else
			t_predSeq{i} = v;
		end
	else
		t_predSeq{i} = '-';
	end
end
phragSeq = [t_predSeq{:}];

bp = phrags(2:end);
end
