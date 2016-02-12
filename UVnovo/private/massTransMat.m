function [massTransitionMat, AAs_tm] = massTransMat(aaCountsFile, paramsAAs)
% MASSTRANSMAT marginalize aa-based transition matrix by unit masses.
% 
% @TODO Rename function & vars. Add documentation.

% This assumes all amino acids and ptms have their own single-char symbol.

if ~exist('paramsAAs','var'), paramsAAs = []; end
AAs = MSaalist(paramsAAs, 'masstype', 'nominal'); % get nominal-mass AA props

% Load file of AA pair and N-|C-term residue counts.

% There's no validation that the data is formatted correctly.
% @TODO bring the code for in silico digest & residue pair counting into repo.
x = importdata(aaCountsFile);
% AA <--> AA transition counts in the forward and backward directions.
tmAA_counts = x.data(1:(end-2), :);
ntermAA_counts = x.data(end-1, :);
ctermAA_counts = x.data(end, :);
% Amino acid symbols. These correspond to x & y dims of tmAA_counts.
aaSyms = strsplit( strtrim( x.textdata{1} ));


[i,ia] = ismember([aaSyms{:}], AAs.aas);
assert(all(i), 'Undefined residue symbols ''%s'' in aa counts file.', ...
	strjoin(aaSyms(~i), ', '))
[aaIntmass,ic] = sort(AAs.aamass(ia)); % sort by ascending mass

aa = [aaSyms{ic}];

% Map between integer mass and symbol(s).
[unimass, uni2full_first, full2uni] = unique(aaIntmass);
uni2full = accumarray(full2uni, 1:numel(aa), [numel(unimass), 1], @(x){x});

t = cellfun(@(x) aa(x), uni2full, 'UniformOutput', false);
mass2symbol = containers.Map(unimass, t);

% Combine isobaric counts in matrix of AA transition counts.
% Frequency of C-terminal residues.
mstart = accumarray(full2uni, ntermAA_counts) / sum(ntermAA_counts);
% Frequency of N-terminal residues.
mend = accumarray(full2uni, ctermAA_counts) / sum(ctermAA_counts);

[xx,yy] = meshgrid( full2uni );
combinedCounts = accumarray({xx(:),yy(:)},tmAA_counts(:));

% Forward prob observing md(i+1) given md(i). [md: mass delta]
mforward  = bsxfun(@rdivide, combinedCounts', sum(combinedCounts,2)');
% Backward prob observing md(i) given md(i+1).
mbackward = bsxfun(@rdivide, combinedCounts,  sum(combinedCounts));

massTransitionMat = struct( ...
	'forward', mforward, ...
	'backward', mbackward, ...
	...'raw', combinedCounts, ...
	'start', mstart, ...
	'end', mend, ...
	'massSyms', aa( uni2full_first ), ... % unique symbol for each mass
	'masses', unimass, ...
	'mass2symbol', mass2symbol);

% Get AAs_tm: a struct of the residues that define valid HMM transistions.
excludeAAs = AAs.aas( ~ismember(AAs.aas, [aaSyms{:}]) );
if any(excludeAAs)
	if isfield(paramsAAs, 'excludeAAs')
		paramsAAs.excludeAAs = [paramsAAs.excludeAAs, excludeAAs];
	else
		paramsAAs.excludeAAs = excludeAAs;
	end
	AAs_tm = MSaalist(paramsAAs, 'masstype', 'nominal');
else
	AAs_tm = AAs;
end





