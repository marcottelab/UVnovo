function params = user_params_NoPar()
% User defined parameters.
%	Any specified here replace the values located in ./UVnovo/default_params().

params.UVnovo.useParallel = false;
% Primary parent mass field name.
params.pre.import_spectra.pmass_field = 'pmass_theorNormmass';


% Amino acids
% @TODO I'm not happy with how this works. Make it better.
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




params.train.ensemble.iters = cell2struct({
	...<nvars>  <ntrees>  <minLeaf>  <oobVarImp>  <plotOOBError>
		0,       8,        1,        'off',        true;    % round 1
		60,      8,        2,        'off',        false;   % round 2
		30,      10,       1,        'on',         true;    % final round
	}, ...
	{  'nvars', 'ntrees', 'minLeaf', 'oobVarImp', 'plotOOBError'}, 2);

