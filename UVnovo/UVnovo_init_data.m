function UVnovo_init_data(spectra_filepath, psm_filepath, paramsIn)
% UVNOVO_INIT_DATA imports, partitions, and saves spectra & psm data.
% 
%	Vars are saved in a serialized matlab format. The is much more efficient in
%	time and space than standard matlab files, but io.loadSer() is required for
%	opening such files. See io.saveSer for more info.
% 


%% Define experiment files & directory structure.

% Prompt for MS2 spectrum file if not provided or accessible.
if ~exist('spectra_file','var') || (exist(spectra_filepath,'file') ~= 2)
	[fn, pn] = uigetfile('*.ms2', 'select ms2 file');
	spectra_filepath = fullfile(pn, fn);
	if ~pn, return, end
end
[pn_ms2, fn_ms2] = fileparts(spectra_filepath);

% Prompt for PSMs file if not provided or accessible.
if ~exist('fn_psms','var') || (exist(psm_filepath,'file') ~= 2)
	[fn, pn] = uigetfile('*.xlsx;*.xls', 'select psm file', pn_ms2);
	psm_filepath = fullfile(pn, fn);
	if ~pn, return, end
end
[~, fn_psm] = fileparts(psm_filepath);

% Construct 'Meta' struct. It will contain paths & parameters for current work.
Meta.paths.spectra = struct( ...
		'name', fn_ms2, ...
		'path', spectra_filepath ...
	);

Meta.paths.psms = struct( ...
		'name', fn_psm, ...
		'path', psm_filepath ...
	);

Meta.paths.basedir = pn_ms2;
Meta.experiment = [Meta.paths.spectra.name, datestr(now,'_yyyymmdd')];

% Create a directory for current experiment.
% @TODO allow user option of choosing path & experiment name.
t_pnExp = fullfile(Meta.paths.basedir, '_experiments', Meta.experiment);
expdir = io.genUniquePath( t_pnExp ); % Give it a unique name.
mkdir(expdir)

Meta.paths.expdir = expdir;


%% Parameters

% Get default parameters and update with any user-defined params.
%	User params are defined in either or both of 'paramsIn' and ./user_params.m,
%	and 'paramsIn' takes precedence. See init_params() for 'paramsIn' format.
t_params = init_params(default_params, user_params);

if ~exist('paramsIn','var'), paramsIn = []; end
Meta.params = init_params(t_params, paramsIn);


%% Import spectra and psm files.

msData = import_spectra(Meta.paths.spectra.path, Meta.params.pre.import_spectra);
psmData = import_psms(Meta.paths.psms.path);


%% Partition and save.

fprintf(1,'Constructing training and test sets and saving to dir:\n\t%s\n',expdir)

% Partition for cross-validation.
nPartitions = Meta.params.pre.crossVal.nPartitions;
partitions = partition_psms(psmData, nPartitions);

% Index between psmData.scans & msData.scans. Assign later to same partitions.
mapi = cross_index(psmData, msData);
for i = 1:nPartitions
	partitions(i).test.msData =  mapi.psmData.msData(partitions(i).test.psmData);
	partitions(i).train.msData = mapi.psmData.msData(partitions(i).train.psmData);
end

% Create and save training and test sets.
%	Initialize variables for saving. We replace 'msData' and 'psmData' with
%	partition-specific data so that the variables are named consistently upon
%	reloading the data.
if ~exist('msData_ORIG', 'var'), msData_ORIG =  msData;  end
if ~exist('psmData_ORIG','var'), psmData_ORIG = psmData; end
msData.scans =  [];
psmData.scans = [];


for n = 1:nPartitions
	% Create and save test partition vars.
	msData.partition.n = n;
	msData.partition.type = 'test';
	msData.scans = msData_ORIG.scans( partitions(n).test.msData );
	
	fn_test = fullfile(expdir, sprintf('p%d_Test.mat', n));
	io.saveSer(fn_test, 'msData', 'Meta')
	
	% Create and save training partition vars.
	msData.partition.type = 'train';
	msData.scans = msData_ORIG.scans( partitions(n).train.msData );
	
	psmData.partition.n = n;
	psmData.partition.type = 'train';
	psmData.scans = psmData_ORIG.scans( partitions(n).train.psmData );
	
	fn_train = fullfile(expdir, sprintf('p%d_Train.mat', n));
	io.saveSer(fn_train, 'msData', 'Meta')
end

fprintf(1,'Done.\n')


