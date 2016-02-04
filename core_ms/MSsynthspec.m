function synth_spec = MSsynthspec(seqIn, AAs, ions, prettyprint)
% SYNTH_SPEC = MSSYNTHSPEC(SEQIN, AAS, IONS, PRETTYPRINT)
% calculate fragment masses for an amino acid sequence
% can include specific mass anywhere in sequence to be included on either
%   preceding residue or as independent from a residue
%   ex: 'PEP[+123.456]TIDE' would put an additional 123.456 on the second P
%       'PEP[123.456]TIDE' would count 123.456 as an independent mass between
%			PEP and TIDE
% 
% If AAs (structure of amino acid masses) is 'nominal' mass type, this function
%	converts any mass annotations to nominal integer masses.
% 
% @BUG @TODO this breaks on sequences with >1 consecutive mass annotations where
%	at least one of those is an 'independent' mass.
%	ex. These to cause an error 'PEP[+12][3]T' and 'PEP[12][3]T'.
%		This one is fine 'PEP[+12][+3]T'.
% @BUG @TODO abc,xyz ions are assumed protonated. Masses should be calculated
%	with addition of mProton, not mH. (ex. y_ion = M + mH2O + mProton)
%	Some other packages do the same, but the mass off by 0.000549 Da.
% @TODO refactor input args into standard framework.
% 
% 


%%% Parse input args
	% amino acids and masses
if ~exist('AAs','var') || isempty(AAs), AAs = MSaalist; end
	% fragment ions to create (upper case for charge 2)
if ~exist('ions','var') || isempty(ions), ions = 'by'; end
	% print a nicely formatted table of sequence ion masses
if ~exist('prettyprint','var') || isempty(prettyprint)
	if nargout == 0
		prettyprint = true;
	else
		prettyprint = false;
	end
end


if isempty(seqIn) || strcmp(seqIn,'x')
    synth_spec = zeros(0,numel(ions));
    return
end

%%% allow for mass modifiers in sequence input itself
[massAdditionStrs, mass_starts] = ...
	regexp(seqIn,'(?<=\[)\+?\-?\d*(\.\d*)?(?=\])','match','start');
strs = regexp(seqIn,'((?<=^|\])[a-zA-Z]*)','match');


t = zeros(size(seqIn));
t(regexp(seqIn,'[a-zA-Z]')) = 1;
t = cumsum(t);
massAdd_locations = t(mass_starts);

massAddPTMinds = cellfun(@(x)x(1)=='+'||x(1)=='-', massAdditionStrs);
if ~isempty(massAdd_locations) && ~massAdd_locations(1) && massAddPTMinds(1)
	warning('MSsynthspec:seqInBadNtermMass', ...
		'Nterm seqIn mass defined as a ptm but not attached to any residue')
	massAddPTMinds(1) = false;
end

sequence = strcat(strs{:});
n_aas = size(sequence,2);
%	number of amino acids in sequence plus independent masses
n_masses = n_aas + nnz(~massAddPTMinds);

%	vector of added mass at each location
if ~isempty(massAdd_locations)
	massAdditions = str2double(massAdditionStrs);
	if strcmp(AAs.masstype, 'nominal')
		massAdditions = round(massAdditions./AAs.m.unit_g);
	end
	massAdds_sitewise = accumarray( massAdd_locations'+1, massAdditions, [n_aas+1,1])';
else
	massAdds_sitewise = zeros(1, n_aas + 1);
end

massAddMask = false(1, n_aas + 1);
massAddMask(massAdd_locations(massAddPTMinds)+1) = true;

aa_masses = zeros(1,n_aas);
for i = 1:n_aas
    aa_masses(i) = AAs.aamass(AAs.aas==sequence(i));
end

t = [ [0, aa_masses] + massAdds_sitewise.*massAddMask; ...
	  massAdds_sitewise.*(1-massAddMask) ];
t(t==0)=[];
masses = t;


n_term = cumsum(masses(1:end));
c_term = cumsum(masses(end:-1:1));
synth_spec = zeros(n_masses, numel(ions));
for i = 1:numel(ions)   % AAs.ncderiv is derivitized [N, C] terminus mass
    switch ions(i)
        case {'n','N'}    % amino terminal aa chain
			% peptide N-term mass has an additional Hydrogen, not included here
            synth_spec(:,i) = n_term + AAs.ncderiv(1);
        case {'o','O'}    % carboxy terminal aa chain
			% peptide C-term mass has an additional OH, not included here
            synth_spec(:,i) = c_term + AAs.ncderiv(2);
        case {'a','A'}    % -27 (+H -CO)
            synth_spec(:,i) = n_term + AAs.m.h - AAs.m.co + AAs.ncderiv(1);
        case {'b','B'}    % +1 (+H)
            synth_spec(:,i) = n_term + AAs.m.h + AAs.ncderiv(1);
        case {'c','C'}    % +18 (+H +NH +H +H)
            synth_spec(:,i) = n_term + 3*AAs.m.h + AAs.m.nh + AAs.ncderiv(1);
        case {'x','X'}    % +45 (+OH +CO)
            synth_spec(:,i) = c_term + AAs.m.oh + AAs.m.co + AAs.ncderiv(2);
        case {'y','Y'}    % +19 (+H2O +H)
            synth_spec(:,i) = c_term + AAs.m.h2o + AAs.m.h + AAs.ncderiv(2);
        case {'z','Z'}    % +2 (+OH -NH)
            synth_spec(:,i) = c_term + AAs.m.oh - AAs.m.nh + AAs.ncderiv(2);
    end
    if double(ions(i))<97 % charge 2+ ions
        synth_spec(:,i) = (synth_spec(:,i) + AAs.m.prot)/2;
    end
end

if prettyprint % print table of residues, masses, and cumulative masses
	ss=[zeros(1,size(synth_spec,2)); synth_spec];
	ionHeaders = cell(1,numel(ions)+1);
	ionHeaders{1} = 'residue';
	for i = 1:numel(ions)
		switch ions(i)
			case {'n','N'}
				ionHeaders{i+1} = 'N-term';
			case {'o','O'}
				ionHeaders{i+1} = 'C-term';
				ss(:,i) = flipud(ss(:,i));
			case {'a','A'}
				ionHeaders{i+1} = 'a-ions';
			case {'b','B'}
				ionHeaders{i+1} = 'b-ions';
			case {'c','C'}
				ionHeaders{i+1} = 'c-ions';
			case {'x','X'}
				ionHeaders{i+1} = 'x-ions';
				ss(:,i) = flipud(ss(:,i));
			case {'y','Y'}
				ionHeaders{i+1} = 'y-ions';
				ss(:,i) = flipud(ss(:,i));
			case {'z','Z'}
				ionHeaders{i+1} = 'z-ions';
				ss(:,i) = flipud(ss(:,i));
		end
		if double(ions(i))<97
			ionHeaders{i+1} = [ionHeaders{i+1} '+2'];
		end
	end
	
	if strcmp(AAs.masstype, 'nominal')
		pre = regexp(massAdditionStrs,'^[+-]','match','once');
		massAdditionStrs = strcat(pre, sprintfc('%g',massAdditions));
	end
	
	z = cellstr(sequence');
	x = cell(n_masses,1);
	xi = 1;
	zi = 1;
	for r = 1:numel(massAdd_locations)
		i = massAdd_locations(r);
		while zi <= i
			x{xi} = z{zi};
			xi = xi+1;
			zi = zi+1;
		end
		if massAddPTMinds(r)
			x{xi-1} = [x{xi-1} '[' massAdditionStrs{r} ']'];
		else
			x{xi} = ['[' massAdditionStrs{r} ']'];
			xi = xi+1;
		end
	end
	x(xi:end) = z(zi:end);
	
	
	maxlen = max(length(ionHeaders{1}), max(cellfun('length',x)));
	frmtSpec =  cell(1,numel(ions)+1);
	frmtSpecHdr =  cell(1,numel(ions)+1);
	frmtSpec{1} = ['\t%-' num2str(maxlen) 's'];
	frmtSpec(2:1+numel(ions)) = {'%-9.3f'};
	frmtSpecHdr{1} = ['\n\t%-' num2str(maxlen) 's'];
	frmtSpecHdr(2:1+numel(ions)) = {'%-9s'};
	
	ntermDiffs = [];
	t = find(ismember(ions,'nabc'),1);
	if t
		ntermDiffs = num2cell([diff(ss(:,t)); 0]);
		ionHeaders{end+1} = 'diff_N';
		frmtSpec{end+1} = '%-6.0f';
		frmtSpecHdr{end+1} = '%-6s';
	end
	ctermDiffs = [];
	t = find(ismember(ions,'oxyz'),1);
	if t
		ctermDiffs = num2cell([0; -diff(ss(:,t))]);
		ionHeaders{end+1} = 'diff_C';
		frmtSpec{end+1} = '%-6.0f';
		frmtSpecHdr{end+1} = '%-6s';
	end
	
	seqMasses = [[x;'-'], num2cell(ss), ntermDiffs, ctermDiffs]';
	seqMassDisp = [sprintf( [strjoin(frmtSpecHdr,'  '),'\n'], ionHeaders{:}),...
		sprintf( [strjoin(frmtSpec,'  '),'\n'], seqMasses{:})];
	
	fprintf(seqMassDisp)
	
	if nargout == 0
		clear synth_spec
	end
end











