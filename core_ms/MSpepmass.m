function seqmass = MSpepmass(sequence, aas, mhPlus, aaOnly)
% MSPEPMASS() calculate peptide mass from sequence(s)
% INPUT ARGS
%   sequence: <str|cellstr> peptide sequences
%       Each can contain additional mass annotations within brackets.
%       Ex: PE[-45.45]TID[.1]E[+55]
%   aas: <struct> from MSaalist, struct of masses for peptide mass calculation
%   mhPlus: true (default) if requesting MH+ mass (addition of a proton)
%   aaOnly: sequence(s) only contain aas symbols. ~10x faster if true.
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

if aaOnly
    if ~iscell(sequence)
		if isempty(sequence)
			seqmass = 0;
		else
			seqmass = sum(AAs.intaa2mass(double(sequence))) + AAs.m.h2o + massProtOrNot + sum(AAs.ncderiv);
		end
	else
		seqmass = zeros(size(sequence));
		intaa2mass = AAs.intaa2mass;
		m_const = AAs.m.h2o + massProtOrNot + sum(AAs.ncderiv); % mass added to all seqs
		for i = 1:numel(sequence)
			s = sequence{i};
			if ~isempty(s)
				seqmass(i) = sum(intaa2mass(double(s))) + m_const;
			end
		end
    end
else % allow PTM annotations
	intaa2mass = AAs.intaa2mass;
	intaa2mass(double('[]')) = 0;
	m_const = AAs.m.h2o + massProtOrNot + sum(AAs.ncderiv);
    if ~iscell(sequence)
		if isempty(sequence)
			seqmass = 0;
		else
			[seqParts, doubMass] = regexp(sequence, ...
				'(?<=\[).*?(?=\])', 'split', 'match');
			seq = [seqParts{:}];
			userdefMass = sum(str2double(doubMass));
			
			seqmass = sum(intaa2mass(double(seq))) + m_const + userdefMass;
		end
	else
		seqmass = zeros(size(sequence));
        for i = 1:numel(sequence)
			s = sequence{i};
			if ~isempty(s)
				[t_seqParts, t_doubMass] = regexp(s, ...
					'(?<=\[).*?(?=\])', 'split', 'match');
				t_seq = [t_seqParts{:}];
				
				%	This is much faster than builtin str2double().
				[~,ia] = max(cellfun('length', t_doubMass));
				t_doubMass{ia}(end+1) = ' ';
				t_userdefMass = sum( sscanf( char(t_doubMass)', '%f'));
				% t_userdefMass = sum(str2double(t_doubMass));
				
				seqmass(i) = sum(intaa2mass(double(t_seq))) + m_const + t_userdefMass;
			end
        end
    end
end

end

