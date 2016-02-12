function params = user_params()
% User defined parameters.
%	Any specified here replace the values located in ./UVnovo/default_params().


% Primary parent mass field name.
params.pre.import_spectra.pmass_field = 'pmass_theorNormmass';


% Amino acids
params.AAs = struct( ...
	...	N/C-terminal derivatization fixed mass
	'ncderiv', [CONSTS.mPTM.AMCA, 0], ...
	...	PTMs: <2 x n> cell of { symbol<char>, mass<real>; ... }
	'ptms', ...
		{{
			'm', CONSTS.mAA.M + CONSTS.mPTM.Oxidation;
			'k', CONSTS.mAA.K + CONSTS.mPTM.Carbamyl;
		}} ....
	);
