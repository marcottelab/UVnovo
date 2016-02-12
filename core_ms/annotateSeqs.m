function seqsAnno = annotateSeqs(seqsIn, ptmsIn, ptmLib)
% annotateSeqs: annotate sequences with PTMs and N-Term addition
% 
% INPUT ARGS
%	seqsIn and ptmsIn are both <cellstr [n x 1]> OR <cell [n x 1] of {cellstr}>.
% 	seqsIn: peptide sequences
% 	ptmsIn: peptide modifications, SEQUEST format
%	ptmLib: <struct> modification name -> mass
% OUTPUT
%	seqsAnno: <cellstr [n x 1]> sequences annotated in accordance with
%		MSpepmass, MSsynthspec, and the like.
% 
% EXAMPLES
%	seqsIn = {'mADmAX', 'PEpPHAW'};
%	ptmsIn = {'N-Term(AMCA); M4(Oxidation)', 'P4(Oxidation)'}
%	-->  seqsAnno = {'M[+215.058243]ADM[+15.994915]AX', 'PEpP[+15.994915]HAW'}
%	% NOTE: The first sequence and ptm spec match in that lowercase symbols and
%		ptms align. This yields a well-formatted annotated sequence. The ptm
%		spec for the second sequence does not align with the lower-case Pro, and
%		this symbol remains lowercase and unmodified in the output.
% 
% @TODO This doesn't yet interpret all Sequest modification syntax
% 
% Inner workings: Each unique modification spec string is turned into a regular
%	expression and replacement expression. It's then applied to the unique set
%	of sequences that it originally paired with.
%	This method is neat but not as efficient as possible.
% 
% !! What should I do about 'x' in sequence & 'X3(L)' in mod list??? Doesn't
%	fit into nice framework below.
% 
% @TODO review this function code.
% @TODO Fix edge case. PTM string assumed to match seq when mods stack.
%	ex. 'kMAD', 'N-Term(AMCA);M1(Oxidation)' still adds +215,+16 to k.
% 
% See also MSAALIST, MSPEPMASS, MSSYNTHSPEC.


if (numel(seqsIn) == 0) || ((numel(seqsIn) == 1) && isempty(seqsIn{1}))
	seqsAnno = seqsIn;
	return
end

if ~exist('ptmLib','var') || ~isstruct('ptmLib')
	% struct of ( ptm_name {str} => monoisotopic_mass {double} )
	ptmLib = CONSTS.mPTM;
end

% Undocumented func sprintfcf allocates each printed string to a new cell.
ptmMass_strs = regexprep(sprintfc('%f',cell2mat(struct2cell(ptmLib)),true),'i$','');

% @TODO containers.Map class is slow. See Yair's book for alternative.
ptms = containers.Map(fieldnames(ptmLib), ptmMass_strs(1:end));

if iscellstr(seqsIn)
	seqsIn = num2cell(seqsIn);
	ptmsIn = num2cell(ptmsIn);
	cellstr_input = true;
else
	cellstr_input = false;
end

assert( isequal(cellfun('prodofsize',seqsIn), cellfun('prodofsize',ptmsIn)))


%%% find unique set of seq:mod pairs
seqs = cat(1,seqsIn{:});
mods = cat(1,ptmsIn{:});
mods(~cellfun(@ischar,mods)) = {''};

[s_uni, ~, ic_s] = unique(seqs);
[m_uni, ~, ic_m] = unique(mods);

[a, ~, ic] = unique([ic_s, ic_m], 'rows');

%%% make each mod declaration into a regular expression
t = regexprep(m_uni,'(?<=^|;)N-Term(?=\()','.1');
n = regexp(t, ...
	'(?:^|;)\s*(?<aa>[\.a-zA-Z](?=\d))?(?<pos>1|\d+)\((?<type>.+?)\)','names');

nModVars = numel(m_uni);
modExpr = cell(nModVars,1);
modExprRep = cell(nModVars,1);

valid_ptmTypes = keys(ptms);

for i = 1:nModVars
	nAnnotations = numel(n{i});
	
	if isempty(m_uni{i})
		% do nothing
		
	elseif numel(n{i}) == 0  || ...
			~all(ismember({n{i}.type},valid_ptmTypes))
		% Invalid PTM string or unknown PTM type.
		modExpr{i} = '.*';
		modExprRep{i} = sprintf('Uninterpretable PTM str ''%s''', m_uni{i});
	
	else
		% PTM string was parseable and each is specified in ptmLib.
		pos = sscanf(sprintf('%s ','0',n{i}.pos),'%f');
		aa = {n{i}.aa};
		ptmMass = values(ptms,{n{i}.type});
		
		nsteps = diff(pos)-1;
		multimod = nsteps == -1;
		if any(multimod)
			nsteps(multimod) = 0;
			aa{multimod} = '';
		end
		t = [num2cell(nsteps)'; aa];
		modExpr{i} = ['(?i)^', sprintf('(.{%d}%s)', t{:})];
		
		t = [num2cell(1:nAnnotations); ptmMass];
		modExprRep{i} = sprintf('$%d[%s]', t{:});
	end
end

nSeqMod = size(a,1);

im = accumarray(a(:,2),1:nSeqMod,[nModVars,1],@(x){x});

t_uniseqAnno = cell(nSeqMod,1);
for i = 1:nModVars
	ism = im{i};
	t_uniseqAnno(ism) = regexprep(s_uni(a(im{i},1)), modExpr(i), modExprRep(i));
end
t_uniseqAnno = regexprep(t_uniseqAnno, '[a-z]\[', '${upper($0)}');

if cellstr_input
	seqsAnno = reshape( t_uniseqAnno(ic), size(seqsIn));
	
else % Put annotated seqs back into original cell layout.
	seqsAnno_expanded = t_uniseqAnno(ic);
	
	%	@TODO extract this stuff to other method
	cell_counts = reshape( cellfun('prodofsize',seqsIn), [], 1);
	ind_cellblock = cumsum([1; cell_counts]);
	seqsAnno = cell(size(seqsIn));
	for i = 1:numel(seqsIn)
		seqsAnno{i} = seqsAnno_expanded(ind_cellblock(i):ind_cellblock(i+1)-1);
	end
end



