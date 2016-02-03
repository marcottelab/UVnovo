function params = user_params_NoPar()
% User defined parameters.
%	Any specified here replace the values located in ./UVnovo/default_params().

params.UVnovo.useParallel = false;
params.pre.import_spectra.pmass_field = 'pmass_theorNormmass';


% Amino acids
% @TODO I'm not happy with how this works. Make it better.
params.AAs = struct( ...
	... disallowed amino acids
	'excludeAAs', 'CK', ...  % cysteine was not alkylated in Ecoli samples :(
	...	N/C-terminal derivatization fixed mass
	'ncderiv', [CONSTS.mPTM.AMCA, 0], ...
	...	PTMs: <2 x n> cell of { symbol<char>, mass<real>; ... }
	'ptms', ...
		{{
			'm', CONSTS.mAA.M + CONSTS.mPTM.Oxidation;
			'k', CONSTS.mAA.K + CONSTS.mPTM.Carbamyl;
		}} ....
	);
