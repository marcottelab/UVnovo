function varargout = UVnovo(exec_mode, varargin)
% UVNOVO initializes the UVnovo environment (parameters, paths, etc) and
%	delegates execution to one of four main workflow functions.
% 
% UVNOVO(EXEC_MODE, VARARGIN)
% 
% 	EXEC_MODE program execution workflow
% 	  'partition'    Import data and create training and test sets.
% 	  'train'        UVnovo training. Make random forest from known spectra.
% 	  'denovo'       Perform de novo sequencing on unknown spectra.
% 	  'benchmark'    Compare de novo results to PSM knowledge.
% 
%	VARARGIN name-value pairs
%	  '-h' or '-help'        Print help and exit.
%	  '-recipe'  <filepath>  text file recipe for automated analysis. @TODO!
%	  '-params'  <filepath>  user-defined parameters function.
%	  '-savedir' <dirpath>   directory in which to save results.
%	  '-ms2'     <filepath>  *.ms2 spectra file.
%	  '-psms'    <filepath>  psms file. See IMPORT_PSMS for valid file formats.
%	  '-train'   <filepath>  training data prepared by UVNOVO PARTITION.
%	  '-test'    <filepath>  test data prepared by UVNOVO PARTITION.
%	  '-rf'      <filepath>  random forest constructed by UVNOVO TRAIN.
%	  '-denovo'  <filepath>  de novo sequencing results from UVNOVO DENOVO.
% 
% @TODO add documentation.
% @TODO make a simple GUI for workflow and file selection.
% 
%	See also UVNOVO_TRAIN, UVNOVO_PARTITION, UVNOVO_DENOVO.

if nargin == 0
	help(which(mfilename('fullpath')))
	return
end

% Add directories associated with UVnovo to the Matlab search path.
updatepath()

% Translate command line arguments into a struct.
args = parseArgs(varargin{:});

% Perform selected execution workflow.
switch exec_mode
	case {'partition', 'partition_data'}
		Meta = partition(args);
		varargout{1} = Meta;
		
	case 'train'
		Meta = train(args);
		varargout{1} = Meta;
		
	case 'denovo'
		Meta = denovo(args);
		varargout{1} = Meta;
		
	case 'benchmark'
		benchmark(args)
		
	otherwise % Run any UVnovo subfunction.
		if ~strcmp(which(exec_mode), which(mfilename))
			error('UVnovo:InvalidExecMode','Arg1 must be a UVnovo subfunction.')
		end
		fh = str2func(exec_mode);
		[varargout{1:nargout}] = fh(varargin{:});
end

end %UVnovo


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%   UVnovo Workflows
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Import data and create training and test sets.
function Meta = partition(argsIn)
	
	% Pertinent commmand line arguments when running UVnovo in partition mode.
	% Any not provided will be queried for, if required for execution.
	argsTemplate = struct(...
		'fn_spectra', '', ...
		'fn_psms',    '', ...
		'fn_params',  uvnovo_rootdir, ...
		'savedir',    ''  ...
		);
	args = initParams(argsTemplate, argsIn);
	
	% Structure 'Meta' centralizes parameters, files, & general information.
	Meta = struct;
	
	% Prompt for any file paths not defined in input args.
	% Struct defines which files are needed & how to prompt for any missing:
	%	file type [str]; path [str|empty]; required [0|1]; uigetfile args.
	filespec = cell2struct( { ...
		'spectra', args.fn_spectra, 1,  {'*.ms2', 'select ms2 file', 'MultiSelect','off'};
		'psms',    args.fn_psms,    1,  {'*.xlsx;*.xls', 'select psm file', 'MultiSelect','off'};
		'params',  args.fn_params,  0,  {{'*params*.m';'*.m'}, '[OPTIONAL] select user parameter file', 'MultiSelect','off'};
		}, {'type', 'path', 'required', 'uigetfile_args'}, 2);
	Meta.paths = getFiles( filespec );
	
	% Get UVnovo parameters & update with user-provided params.
	Meta.params = getParams(Meta.paths.params.path);
	
	% Initialize an 'experiment' name & working dir for current workflow.
	default_expname = [Meta.paths.spectra.name, datestr(now,'_yyyymmdd')];
	Meta = initExperiment(Meta, default_expname, args.savedir);
	
	% Import spectra and psm data, partition, and save.
	fprintf(1, 'Constructing training and test datasets and saving at:\n\t%s\n', ...
		Meta.paths.exp.path)
	
	if ~isfield(args,'dryrun') || ~args.dryrun
		Meta.partition = UVnovo_partition(Meta);
	end
	
	fprintf(1,'Done.\n')
end


%%% UVnovo model training. Construct random forests from spectra & psms.
function varargout = train(argsIn)
	
	argsTemplate = struct(...
		'fn_train',   '', ...
		'fn_params',  uvnovo_rootdir, ...
		'savedir',    ''  ...
		);
	args = initParams(argsTemplate, argsIn);
	
	% Structure 'Meta' centralizes parameters, files, & general information.
	Meta = struct;
	
	% Prompt for any file paths not defined in input args.
	% @TODO also allow selection of spectra & psm files instead of training set.
	filespec = cell2struct( { ...
		'train',   args.fn_train,   1,  {{'*Train*.mat';'*.mat'}, 'select training dataset', 'MultiSelect','off'};
		'params',  args.fn_params,  0,  {{'*params*.m';'*.m'}, '[OPTIONAL] select user parameter file', 'MultiSelect','off'};
		}, {'type', 'path', 'required', 'uigetfile_args'}, 2);
	Meta.paths = getFiles( filespec );
	
	% % % Load Meta from training dataset.
	% % s = io.loadSer(Meta.paths.train.path);
	
	% Get UVnovo parameters & update with user-provided params.
	Meta.params = getParams(Meta.paths.params.path);
	
% 	% Update 'Meta' with experiment name & working/output directory.
% 	t_expname = []; % @TODO populate with something.
% 	if ~isempty(args.savedir)
% 		t_expdir = args.savedir;
% 	else
% 		t_expdir = fileparts(Meta.paths.train.path);
% 	end
% 	Meta = initExperiment(Meta, t_expname, t_expdir);
	
	% Open matlabpool if using parallel computation & it's not already open.
	% A just-opened pool closes when 'closepool_onCleanup' is destroyed, i.e.
	% when this function returns or there's an error somewhere.
	closepool_onCleanup = openParallel(Meta); %#ok<NASGU>
	
	% Build models from training data.
	fprintf(1, 'Training UVnovo models ... \n')
	
	if ~isfield(args,'dryrun') || ~args.dryrun
		UVnovo_train(Meta)
	end
	
	varargout = {Meta};
	fprintf(1,'All done.\n')
end


%% Perform de novo sequencing on unknown spectra
function varargout = denovo(argsIn)
	
	argsTemplate = struct(...
		'fn_spectra', '', ...
		'fn_rf',      '', ...
		'fn_AAmodel', fullfile(uvnovo_rootdir,'UVnovo','models','aaPairCounts_Ecoli_K12.txt'), ...
		'fn_params',  uvnovo_rootdir, ...
		'savedir',    ''  ...
		);
	args = initParams(argsTemplate, argsIn);
	
	% Structure 'Meta' centralizes parameters, files, & general information.
	Meta = struct;
	
	% Prompt for any file paths not defined in input args.
	% @TODO allow selection of spectra file instead of test set.
	filespec = cell2struct( { ...
		'spectra', args.fn_spectra,1, {{'*.ms2';'*Test*.mat';'*.mat'}, 'select *.ms2 OR test set file', 'MultiSelect','off'}; 
		'rf',      args.fn_rf,     1, {{'*RF*.mat';'*.mat'}, 'select trained fragmentation model', 'MultiSelect','off'}; 
		'aamodel', args.fn_AAmodel,1, {{'aaPairCounts*.mat';'*.mat'}, 'select trained fragmentation model', 'MultiSelect','off'}; 
		'params',  args.fn_params, 0, {{'*params*.m';'*.m'}, '[OPTIONAL] select user parameter file', 'MultiSelect','off'}; 
		}, {'type', 'path', 'required', 'uigetfile_args'}, 2);
	
	% Differentiate between test.mat and *.ms2 spectra files.
	% @TODO this should be elsewhere. See getFiles_refactoring.
	t_paths = getFiles( filespec );
	[~,~,ext] = fileparts(t_paths.spectra.path);
	if strcmpi(ext,'.mat')
		t_paths.test = t_paths.spectra;
		t_paths = rmfield(t_paths, 'spectra');
	end
	Meta.paths = t_paths;
	
	% Get UVnovo parameters & update with user-provided params.
	Meta.params = getParams(Meta.paths.params.path);
	
% 	% Update 'Meta' with experiment name & working/output directory.
% 	t_expname = []; % @TODO populate with something.
% 	if ~isempty(args.savedir)
% 		t_expdir = args.savedir;
% 	else
% 		t_expdir = fileparts(Meta.paths.test.path);
% 	end
% 	Meta = initExperiment(Meta, t_expname, t_expdir);
	
	% Open matlabpool if using parallel computation & it's not already open.
	closepool_onCleanup = openParallel(Meta); %#ok<NASGU>
	
	% De novo sequencing!
	fprintf(1, 'De novo sequencing ... \n')
	
	if ~isfield(args,'dryrun') || ~args.dryrun
		UVnovo_denovo(Meta)
	end
	
	varargout = {Meta};
	fprintf(1,'Done.\n')
end


%% Compare de novo results to PSM knowledge
function varargout = benchmark(argsIn)
	
	argsTemplate = struct(...
		'fn_denovo',  '', ...
		'fn_psms',    '', ...
		'fn_params',  uvnovo_rootdir, ...
		'savedir',    ''  ...
		);
	args = initParams(argsTemplate, argsIn);
	
	% Structure 'Meta' centralizes parameters, files, & general information.
	Meta = struct;
	
	% Prompt for any file paths not defined in input args.
	filespec = cell2struct( { ...
		'denovo',     args.fn_denovo, 1, {{'*Predictions*.mat';'*.mat'}, 'select de novo results', 'MultiSelect','off'};
		'psms',       args.fn_psms,   1, {'*.xlsx;*.xls', 'select psm file', 'MultiSelect','off'};
		'params',     args.fn_params, 0, {{'*params*.m';'*.m'}, '[OPTIONAL] select user parameter file', 'MultiSelect','off'};
		}, {'type', 'path', 'required', 'uigetfile_args'}, 2);
	Meta.paths = getFiles( filespec );
	
	% Get UVnovo parameters & update with user-provided params.
	Meta.params = getParams(Meta.paths.params.path);
	
	% Update 'Meta' with experiment name & working/output directory.
	t_expname = []; % @TODO populate with something.
	if ~isempty(args.savedir)
		t_expdir = args.savedir;
	else
		t_expdir = fileparts(Meta.paths.denovo.path);
	end
	Meta = initExperiment(Meta, t_expname, t_expdir);
	
	% Benchmarking
	fprintf(1, 'Benchmarking UVnovo predictions ... \n')
	
	UVnovo_benchmark(Meta)
	
	varargout = {Meta};
	fprintf(1,'Done.\n')
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%   UTILITY FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function uvnovo_root = uvnovo_rootdir()
	uvnovo_root = fileparts(mfilename('fullpath'));
end

function updatepath()
	% Add directories associated with UVnovo to the Matlab search path.
	append_relpaths = {
		'./'
		'./core_ms'
		'./core_ms/fileIO'
		'./UVnovo'
		'./UVnovo/models'
		'./UVnovo/plotting'
		'./UVnovo/to_deprecate'
		'./utils'
		};
	append_paths = fullfile(uvnovo_rootdir, append_relpaths);
	addpath(append_paths{:})
	addpath(genpath(fullfile(uvnovo_rootdir,'dependencies')))
end

function args = parseArgs(varargin)
	args = struct;
	
	% These flags specify files that can be used during execution.
	%	Any file names passed this way are return in args.(varname).
	fileflag_varname = {
	... <flag>     <varname>
		'-recipe', 'fn_recipe'
		'-params', 'fn_params'
		'-ms2',    'fn_spectra'
		'-psms',   'fn_psms'
		'-train',  'fn_train'
		'-test',   'fn_spectra'
		'-rf',     'fn_rf'
		'-denovo', 'fn_denovo'
		};
	
	fileflags = fileflag_varname(:,1);
	flag2var = containers.Map(fileflags, fileflag_varname(:,2));
	
	indFlags = find(~cellfun('isempty', regexp(varargin,'^-')));
	
	% Values beyond the first for each flag are ignored.
	indVals = min(indFlags+1, [indFlags(2:end)-1, numel(varargin)]);
	
	nFlags = numel(indFlags);
	unkFlags = {};
	for n = 1:nFlags
		flag = varargin{indFlags(n)};
		val = varargin{indVals(n)};
		
		switch lower(flag)
			
			% Print help and exit.
			case {'-h','-help'}
				help(which(mfilename('fullpath')))
				error('UVnovo:NoFileSelected', 'No file selected.')
			
			% Files
			case fileflags
				args.(flag2var(flag)) = val;
			
			% Working/save directory
			case {'-savedir'}
				args.savedir = val;
			
			% Dry run. Exit before any intensive processing.
			case {'-dryrun'}
				args.dryrun = true;
			
			% Unknown
			otherwise
				unkFlags = [unkFlags flag]; %#ok<AGROW>
		end
	end
	if ~isempty(unkFlags)
		error('UVnovo:unkArg','Cannot parse argument flag(s): %s', ...
			strjoin(unkFlags,' '))
	end
end

function files = getFiles(filespec)
	% Prompt for each file in 'filespec' if it is not already specified.
	% 
	% filespec <struct>
	%	.type: <str (valid field name)> type of file, ex. 'spectra', 'psms', ...
	%	.path: <str> filepath or dir. User prompt if not valid file. When dir is
	%		provided, uigetfile prompt opens in that location.
	%	.require: <logical> File required? If not, returns empty string when no
	%		file is selected.
	%	.uigetfile_args: <cell> arguments for uigetfile() if path not defined.
	% 
	% @TODO What if user allowed to select multiple files?
	% @TODO improve input file handling when there's overlap between file
	%	purposes. Ex: [train | ms2/psms], [test|ms2]. See getFiles_refactoring.
	
	files = struct;
	nFiles = numel(filespec);
	
	pwd_pre = pwd;
	z = onCleanup(@() cd(pwd_pre));
	for n = 1:nFiles
		pathIn = filespec(n).path;
		
		if ~isempty(pathIn) && ( exist(pathIn,'file') == 2 )
			% File exists.
			filepath = pathIn;
			
		else% Prompt for file.
			
			% Change uigetfile start dir, using longest existing path from
			% provided string. This will have no effect if pathIn is empty.
			cd( io.validPath( pathIn ))
			
			[fn_ext, pn] = uigetfile( filespec(n).uigetfile_args{:} );
			
			if isequal(fn_ext, 0) % No file selected.
				if filespec(n).required
					error('UVnovo:NoFileSelected', ...
						'Required %s file not provided.', filespec(n).type)
				end
				filepath = '';
			else
				filepath = fullfile(pn, fn_ext);
			end
		end
		
		[pn, fn] = fileparts(filepath);
		files.(filespec(n).type) = struct( ...
			'name', fn, ...
			'path', filepath ...
			);
		if ~isempty(pn)
			cd(pn) % Change dir so subsequent uigetfile's start in relevant dir.
		end
	end
end

function params = getParams(varargin)
	% GETPARAMS Get default parameters and update with any user-defined params.
	%	User params are defined in user_params.m and/or passed in varargin.
	% 
	% VARARGIN: [user_params filepath (optional)] [additional user parameters]
	% PARAMS: <struct> parameters for UVnovo execution.
	% 
	%	See also USER_PARAMS, DEFAULT_PARAMS.
	
	%	@TODO put user params into simple format, not a function.
	%		(Refer to psmParse.m and Yair's propertiesGUI.m.)
	
	
	% Revert path on error of function return.
	pwd_pre = pwd;
	z = onCleanup(@() cd(pwd_pre));
	
	% Load user params from provided file or from UVnovo user_params().
	if nargin && exist(varargin{1}, 'file') && ...
			exist(fileparts(varargin{1}), 'dir') % Verify it's not just the name of a random mfile on path.
		% @TODO add more error checking & protections found in RUN.
		[pn, fn] = fileparts(varargin{1});
		cd(pn)
		paramsUser = feval(fn);
		
		varargin(1) = [];
	else
		% Make sure we're not using an unexpected 'user_params' file.
		cd(uvnovo_rootdir)
		rundir = fileparts(which('user_params'));
		if ~strcmp(rundir, uvnovo_rootdir);
			fprintf(1, 'user_params.m is not in expected location:\n\t%s\n', rundir)
			d = input('Use this file (1)? Use default params (2)? Or exit (0)?: ','s');
			switch d
				case '1'
					paramsUser = user_params;
				case '2'
					paramsUser = [];
				otherwise
					error('UVnovo:NoFileSelected', 'No file selected.')
			end
		else
			paramsUser = user_params;
		end
		
		if nargin && isempty(varargin{1})
			% In case an empty path is passed for user_params.
			% The need for this if block is kind of tedious & ugly.
			varargin(1) = [];
		end
	end
	
	% Load default parameters.
	paramsDef = default_params;
	
	% Update defaults with any differences from user params and varargin.
	params = initParams( initParams(paramsDef, paramsUser), varargin);
end

function Meta = initExperiment(Meta, default_expname, default_savedir)
	
	if ~exist('Meta','var')
		Meta = struct;
	end
	if ~exist('default_expname','var') || isempty(default_expname)
		if isfield(Meta, 'experiment')
			default_expname = Meta.experiment;
		else
			default_expname = datestr(now,'uvnovo_exp_yyyymmdd');
		end
	end
	if ~exist('default_savedir','var')
		try
			default_savedir = Meta.paths.exp.path;
		catch
			default_savedir = uvnovo_rootdir;
		end
	end
	
	% Get experiment name.
	% @TODO validate user-defined name.
	prompt = 'Enter experiment name.';
	expname = inputdlg(prompt, 'UVnovo experiment', ...
		[1, max(30, 3 + length(default_expname))], {default_expname}, 'on');
	
	% Get experiment working directory.
	% @TODO uigetdir is clunky. I have something better than this somewhere.
	expdir = uigetdir(default_savedir, 'Choose output directory for experiment.');
	
	if isequal(expdir, 0)
		error('UVnovo:NoOutputDirChosen', 'No output directory selected.')
	end
	
	Meta.experiment = expname;
	Meta.paths.exp.name = '';
	Meta.paths.exp.path = expdir;
	
	% OLD stuff
	% 	Meta.paths.basedir = pn_ms2;
	% 	Meta.experiment = [Meta.paths.spectra.name, datestr(now,'_yyyymmdd')];
	% 	
	% 	% Create a directory for current experiment.
	% 	t_pnExp = fullfile(Meta.paths.basedir, '_experiments', Meta.experiment);
	% 	expdir = io.genUniquePath( t_pnExp ); % Give it a unique name.
	% 	mkdir(expdir)
end

function closepool_onCleanup = openParallel(Meta)
	
	pool_opened = false;
	if Meta.params.UVnovo.useParallel
		if isempty(ver('distcomp')) % Is parallel computing toolbox available?
			warning('UVnovo:openParallel:DistcompNotAvailable', ...
				'Parallel computing toolbox is not available. UVnovo will run in a single process.')
		
		elseif ~matlabpool('size')
			a = Meta.params.UVnovo.matlabpool_varargs;
			if iscell(a)
				matlabpool(a{:})
			else
				matlabpool(a)
			end
			pool_opened = true;
		end
	end
	
	if nargout
		if pool_opened
			% Close parallel pool when 'closepool_onCleanup' is destroyed.
			closepool_onCleanup = onCleanup(@()matlabpool('close'));
		else
			closepool_onCleanup = 0;
		end
	end
	
end

%%

