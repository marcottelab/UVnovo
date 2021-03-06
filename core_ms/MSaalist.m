function AAs = MSaalist(varargin)
% Initialize amino acids, their masses, and other needed masses.
% 
% Masses are by default exact monoisotopic, measured in Da.
%	When parameter 'masstype' is set to 'nominal', this returns nominal (integer
%	count of protons + neutrons) values. PTM and ncderiv nominal masses are
%	estimated as round(Dalton_mass/CONSTS.unit_g).
% 
% OUTPUT
%	AAs <struct>
%		.masstype
%		.aas
%		.aamass
%		.m
%		.ncderiv
%		.intaa2mass   % <array> Index into this array with peptide string or AA
%						to retrieve residue masses.
% 
% @TODO Validate some of the parameters.
% @TODO tweak exact AA masses so isobaric pairs (ex. N/GG, FS/AY) are equal.
% 	ex. diff(sum(AAs.intaa2mass(['GG';'N ']), 2)) % This should be 0.
% @TODO rewrite this and make it better. -- In progress. See MSaalist_refactor.
%	Mainly, the ptm definitions could be more flexible & powerful.
% 
% See also MSPEPMASS, MSSYNTHSPEC, ANNOTATESEQS.

%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%

paramsDef = struct( ...
	'masstype', 'mono', ... 'mono'|'monoisotopic' or 'nominal'
	'ptms', ... <2 x n> cell of {symbol<char>, mass<real>; ...}
		{{
			'm', CONSTS.mAA.M + CONSTS.mPTM.Oxidation;
		}}, ...
	'ncderiv', [0 0], ... N- or C-terminal fixed masses (derivitization)
	'excludeAAs', '' ... <str> exclude each symbol from output
	);

params = initParams(paramsDef, varargin{:});
assert(numel(params) == 1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%

%	amino acid one-letter abbreviations
t = fieldnames(CONSTS.mAA);
aas = [t{:}];
%	amino acid masses
aamass = struct2array(CONSTS.mAA)';

m.prot = CONSTS.mProton;
m.h = CONSTS.mH;
m.h2o = CONSTS.mH2O;
m.co = CONSTS.mCO;
m.nh = CONSTS.mNH;
m.oh = CONSTS.mOH;
% m.unit_g is assigned below.

ptms = params.ptms;
if ~isempty(ptms)
	for i = 1:size(ptms,1)
		residueSymbol = ptms{i,1};
		residueMass = ptms{i,2};
		aaloc = find(aas == residueSymbol);
		if isempty(aaloc)
			aas(end + 1) = residueSymbol;
			aaloc = numel(aas);
		end
		aamass(aaloc) = residueMass;
	end
end

%	Remove specific residues from output.
indExclude = ismember(aas, params.excludeAAs);
aas(indExclude) = [];
aamass(indExclude) = [];

%	Sort aas, aamass by ascending symbol then mass.
[~, iaa] = sort(aas);
[~, iaam] = sort(aamass(iaa));
aas = aas(iaa(iaam));
aamass = aamass(iaa(iaam));

switch params.masstype
	case {'mono','monoisotopic'}
		ncderiv = params.ncderiv;
		AAs.masstype = 'monoisotopic';
		AAs.aas = aas;
		AAs.aamass = aamass;
		AAs.m = m;
		AAs.ncderiv = ncderiv;
		
	case 'nominal'
		% convert to 'nominal' mass units
		AAs.masstype = 'nominal';
		AAs.aas = aas;
		AAs.aamass = round(aamass./CONSTS.unit_g);
		AAs.m = structfun(@(x)round(x./CONSTS.unit_g), m, 'UniformOutput',false);
		AAs.ncderiv = round(params.ncderiv./CONSTS.unit_g);
		
	otherwise
		error('MSaalist:UnkMassType','Unknown mass type ''%s''',params.masstype)
end

% Mean Da per nucleon. Used for nominal mass conversions. See CONSTS.
AAs.m.unit_g = CONSTS.unit_g;

aakeyvals = nan(127,1);
aakeyvals(' ') = 0;
aakeyvals(aas) = AAs.aamass;
AAs.intaa2mass = aakeyvals;


