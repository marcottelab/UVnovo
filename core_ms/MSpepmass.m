function seqmass = MSpepmass(sequences, aas, mhPlus, aaOnly)
% MSPEPMASS() calculate peptide mass from sequence(s)
% INPUT ARGS
%   sequences: <str|cellstr> peptide sequences
%       Each can contain additional mass annotations within brackets.
%       Ex: PE[-45.45]TID[.1]E[+55]
%   aas: <struct> from MSaalist, struct of masses for peptide mass calculation
%   mhPlus: true (default) if requesting MH+ mass (addition of a proton)
%   aaOnly: sequence(s) only contain aas symbols. A little faster if true.
%       When true, does not allow masses in brackets (used to signify PTMs).
%       ex. if sequence == PEPTM[+15.994915]IDE, add mass of the oxidation.
% OUTPUT
%   seqmass: MH+ mass(es)
%   
%   Will return NaN if error in sequence or unable to calculate!!
% 
%   SEQMASS may vary by eps(mass) depending on order of residues.
%       ex. QQEFVDLIR vs. QQEVFDLIR
% 
%	Just for anonymous fun, here's a one-liner that even allows for the
%	bracketed masses ('seqs' should be a cellstr):
%	f_pepmass = @(seqs,aas) sum( aas.m.h2o + aas.m.prot + sum(aas.ncderiv)) + cellfun(@(x)sum(str2double(x)), regexp(regexprep(seqs,{'(^\[|\]$)','[A-Z]+'},{'','${num2str(sum(aas.intaa2mass($0)),''%.6f'')}'}), '(\]\[|[\[\]])', 'split'));
% 

if ~exist('mhPlus','var')||isempty(mhPlus), mhPlus = true; end
if ~exist('aaOnly','var')||isempty(aaOnly), aaOnly = false; end

persistent AA_default   % list of amino acids and masses
if exist('aas','var') && ~isempty(aas), AAs = aas;
else
    if isempty(AA_default)
        AA_default = MSaalist;
    end
    AAs = AA_default;
end

if mhPlus
    massProtOrNot = AAs.m.prot;  % peptide MH+ mass
else
    massProtOrNot = 0;           % peptide M mass
end
m_const = AAs.m.h2o + massProtOrNot + sum(AAs.ncderiv); % mass added to all seqs

% mass indexed to its aa:
aa2mass = AAs.intaa2mass;

if ~iscell(sequences)
	sequences = {sequences};
end
seqmass = zeros(size(sequences));

if ~aaOnly
	% Allow bracketed PTM annotations.
	[seqParts, bracMasses] = regexp(sequences,'(?<=\[).*?(?=\])', 'split', 'match');
	
	if all( cellfun('isempty',bracMasses) )
		% Then there aren't any (valid) mass annotations. Use simple calc below.
		aaOnly = true;
		
	else
		aa2mass(double('[]')) = 0;
		nominalmass = strcmp(AAs.masstype, 'nominal');
		unit_g = AAs.m.unit_g;
		
		for i = 1:numel(sequences)
			s = sequences{i};
			if ~isempty(s)
				t_brackMassStrs = bracMasses{i};
				
				if isempty(t_brackMassStrs)
					t_brackMass = 0;
					
				else
					% This conversion is much faster than builtin str2double().
					[~,ia] = max(cellfun('length', t_brackMassStrs));
					t_brackMassStrs{ia}(end+1) = ' ';
					t_brackMass = sum( sscanf( char(t_brackMassStrs)', '%f'));
					% SLOWER: t_brackMass = sum(str2double(t_brackMassStrs));
					
					if nominalmass
						t_brackMass = round(t_brackMass./unit_g);
					end
				end
				
				seqmass(i) = sum(aa2mass(double([seqParts{i}{:}]))) + ...
					m_const + t_brackMass;
			end
		end
	end
end

if aaOnly
	% No bracketed PTM annotations.
	for i = 1:numel(sequences)
		s = sequences{i};
		if ~isempty(s)
			seqmass(i) = sum(aa2mass(double(s))) + m_const;
		end
	end
end

end

%%% Note to self--
% Matlab will convert chars to numbers when using them to index
% into and array, but it's a little slower than explicitly
% converting to double:
%	aa2mass(double(seq))   % faster
%	aa2mass(seq)           % slower

