function varargout = cleanws(action, varargin)
%CLEANWS  Manage variables in workspace and keep it tidy.
%	The basic functionality allows a user to specify 'protected' vars (in the
%	base workspace) and selectively clear unprotected vars at future times. This
%	only protects vars from being cleared by this function. It does not affect
%	variable access or modification by other means. Additional capabilities
%	include formated display of protected or unprotected vars, finer control
%	over variable inclusion/exclusion/deletion, and cleanws operation.
% 
%   CLEANWS ON saves names of all vars in workspace into persistent state.
%   CLEANWS ON VAR1 VAR2 ... performs as a CLEANWS OFF VARARGIN, CLEANWS ON.
%		If no saved var state, clears all but varargin and stores those names.
%		If var in varargin in not allocated, that name is NOT retained!**
%			@TODO **inconsistent behavior versus CLEANWS APPENDVARS. Rectify.
% 
%   CLEANWS OFF removes all variables in workspace created after CLEANWS ON
%   CLEANWS OFF VAR1 VAR2 ... removes all new variable not among those specified
% 
%   CLEANWS APPENDVARS VARNAME1 VARNAME2 ... add var names to protected list.
%		This will protect a var name even if it's not present in workspace.
%		Subsequent call to CLEANWS ON ... will remove any protected var names 
%			that are not present in workspace.
% 
%	CLEANWS VARS print list of workspace variables to keep
%	CLEANWS VARSS long form of VARS, includes fields: size, bytes, class, etc.
%		This calls builtin 'whos' to produce output.
%	CLEANWS UNPROTECTEDVARS print list of workspace variables not protected
%	CLEANWS UNPROTECTEDVARSS long form of UNPROTECTEDVARS, as with VARS/VARSS.
% 
%	CLEANWS TOGGLE turns cleanws functionality on or off. When disabled, calls
%	have no effect. Persistent states initialized before toggling off remain
%	constant until function is reenabled. This is useful for disabling calls
%	without commenting out every CLEANWS call in executed scripts.
%	CLEANWS('TOGGLE',TRUE) set cleanws to active. Enable future commands.
%	CLEANWS('TOGGLE',FALSE) cleanws inactive. Future commands have no effect.
% 
%	CLEANWS WARNINGS turn cleanws warnings on or off.
%	CLEANWS('WARNINGS',TRUE) set cleanws to print warnings.
%	CLEANWS('WARNINGS',FALSE) cleanws does not print warnings.
% 
%	CLEANWS RESET re-initialize to fresh cleanws state.
% 
%   s = CLEANWS('ON') returns list of protected variables
%   s = CLEANWS('APPENDVARS') returns list of protected variables
%   s = CLEANWS('OFF','var1','var2',...) returns list of cleared variables
% 
%	CLEANWS BASE ...  allows other functions to call cleanws. This overrides 
%	internal protections and should not be abused. As an example usage case,
%		this allows binding hotkeys to cleanws using EditorMacro.
%		(http://www.mathworks.com/matlabcentral/fileexchange/24615)
%			ex: bind key combination Alt+Ctrl+i to CLEANWS ON
%				EditorMacro('Alt ctrl i', @()cleanws('base','on'),'run');
%			more:
%				EditorMacro('Alt ctrl c', @()cleanws('base','off'),'run');
%				EditorMacro('Alt ctrl a', @()cleanws('base','varss'),'run');
% 
%	After initialization, CLEANWS stores a list of protected vars and a few
%	other settings in persistent vars local to this function. These persistent
%	vars are destroyed with CLEANWS RESET, whenever the function is modified, or
%	upon calling CLEAR ALL, CLEAR FUN, or CLEAR FUNCTIONS.
% 
%	To enable tab-completion of cleanws arguments when calling this function,
%	edit '<matlab_path>/toolbox/local/TC.xml' to include:
%   <binding name="cleanws" ctype="VAR">
%     <arg argn="1" ctype="VAR" value="on off appendVars vars varss unprotectedVars unprotectedVarss reset toggle warnings"/>
%   </binding>
% 
% @TODO make it take regular expressions and/or wildcards.
%	!! clearvars would easily allow for wildcards & simplify other stuff, too !!
% @TODO implement a priority system. Associate each protected var with a
%	priority value. During clear, have option of clearing only those below
%	certain threshold.
%	ex. pass an int when protecting vars or default to <3>
%		'cleanws off 2' would clear everything 2 & lower
%		'cleanws off' would clear everything 4 & lower
% 
% See also CLEAR, CLEARVARS, WHOS, WHOSS.

%%

if nargout, varargout = {[]}; end
nargin_ = nargin;
dbs=dbstack;
if nargin_ > 0 && strcmpi('base',action)
	nargin_ = nargin_-1;
	try
		action = varargin{1};
		varargin(1) = [];
	end
elseif numel(dbs)>1 && ~all(strcmp(dbs(1).name,{dbs.name}))
	% then this was called from another function and shouldn't act.
	return
end

actions={...<name>         <toggle-able>  <requires init>
			'on',               true,          false;
			'off',              true,          true;
			'vars',             true,          true;
			'varss',            true,          true;
			'unprotectedvars',  true,          true;
			'unprotectedvarss', true,          true;
			'appendVars',       true,          false;
			'reset',            false,         false;
			'toggle',           false,         false;
			'warnings'          false,         false;
	};

validArgs = actions(:,1)';
if nargin_ == 0 || ~ischar(action) || ~any(strcmpi(validArgs,action));
	fprintf(1,'Usage: Cleanws [\b<%s>  <args>]\b\n',strjoin(validArgs,', '))
	return
end

nPrintCols = 3;

persistent initialized varlist state active warnings
if isempty(initialized)
	initInternalState
end
% Lock this function in memory --> clear cannot destroy its persistent vars.
% mlock  % !!!	problem if cleanws is called in a function - saves scope of that
%				fnc, not local	!!!

if ~active && any(strcmpi(action, actions([actions{:,2}]',1) ))
	if warnings
		warning('cleanws:inactive',...
			['CLEANWS is inactive. No action taken. '...
			'Use ''cleanws toggle'' to activate.'])
	end
	return
end

if isempty(state) && any(strcmpi(action, actions([actions{:,3}]',1) ))
	if warnings
		warning('cleanws:emptyState',...
			'CLEANWS is not initialized. No action taken.')
	end
	return
end

switch lower(action)
	case 'reset' % reset function state
		initInternalState
		
	case 'toggle' % enable/disable this function, preserving current state
		if ~isempty(varargin)
			active = setTF(varargin{1}, active, action);
		else
			active = ~active;
			de = {'de',''};
			fprintf(1,'cleanws %sactivated\n', de{active+1})
		end
		
	case 'warnings' % turn warnings on or off if they're annoying
		if ~isempty(varargin)
			warnings = setTF(varargin{1}, warnings, action);
		else
			warnings = ~warnings;
			onoff = {'off','on'};
			fprintf(1,'cleanws warnings %s\n', onoff{warnings+1})
		end
		
	case {'vars', 'varss'}
		if ~nargout % print protected variables
			if ~isempty(varlist)
				[varsPresent, varsAbsent] = varsInBase(varlist);
				if ~isempty(varsPresent)
					fprintf(1,'[\b%s]\b\n', 'protected vars')
					if strcmpi(action,'varss')
						% verbose output
						printVerbose(varsPresent)
					else % list variable names only
						printTable(varsPresent, nPrintCols)
					end
				end
				if ~isempty(varsAbsent)
					fprintf(1,'[\b%s]\b\n', 'vars not in workspace')
					printTable(varsAbsent, nPrintCols)
				end
			else
				fprintf(1,'  [\b%s]\b\n', 'no protected vars')
			end
		else % return protected variables
			varargout = {varlist};
		end
		
	case 'on' % set protected variables, clear unprotected
		state = 1;
		if ~isempty(varargin) % call CLEANWS OFF to clear unspecified vars
			evalin('base', [ 'cleanws off ' strjoin(varargin,' ')] )
		end
		varlist = evalin('base','who')';
		%%% alternate method to accomplish same
		% % vvalid = cellfun(@isvarname,varargin);
		% % varlist = unique([varlist, varargin(vvalid)]);
		% % %	remove vars from varlist if they're not present in workspace
		% % %	clear unprotected
		if nargout, varargout = {varlist}; end
		
	case 'off' % clear unprotected variables
		protectedVars = cat(2, varlist, varargin);
		clearThese = getUnprotectedVars(protectedVars);
		% @TODO Why did I loop rather than construct a single 'clear' line?
		for i = 1:numel(clearThese)
			evalin('base', ['clear ' clearThese{i}]);
		end
		if nargout, varargout = {clearThese}; end
		
	case 'appendvars' % add to protected variables, no clearing of unprotected
		state = 1;
		vvalid = cellfun(@isvarname,varargin);
		varlist = unique([varlist, varargin(vvalid)]);
		if nargout, varargout = {varlist}; end
		
	case {'unprotectedvars', 'unprotectedvarss'}
		protectedVars = cat(2, varlist, varargin);
		unprotectedVars = getUnprotectedVars(protectedVars);
		if ~nargout % print unprotected variables
			if ~isempty(unprotectedVars)
				fprintf(1,'[\b%s]\b\n', 'unprotected vars')
				if strcmpi(action,'unprotectedvarss')
					% verbose output
					printVerbose(unprotectedVars)
				else % list variable names only
					printTable(unprotectedVars, nPrintCols)
				end
			else
				fprintf(1,'\t[\bNone]\b\n')
			end
		else % return unprotected variables
			varargout = {unprotectedVars};
		end
end

%% %%% ------ support functions ------ %%% %%

	function initInternalState
		active = true;
		varlist = {};
		state = [];
		warnings = true;
		initialized = true;
	end
end

function val = setTF(arg, val, action)
	% Check user input argument for true/false vals.
	%	If ARG is valid, set VAL to true/false depending on user input.
	%	If not, return original value of VAL & print provided warning message.
	switch arg
		case {'1',1,'y','yes','true', 'True', 'TRUE', true, 'on' }
			val = true;
		case {'0',0,'n','no', 'false','False','FALSE',false,'off'}
			val = false;
		otherwise
			warning('cleanws:unknownArg', ...
				'CLEANWS unknown arg following ''%s''', action)
	end
end

function [varsPresent, varsAbsent] = varsInBase(varlist)
	i = ismember(varlist, evalin('base','who'));
	varsPresent = varlist(i);
	varsAbsent = varlist(~i);
end

function printTable(varlist, ncols) % extracted to io.printTable
	maxl = max(cellfun(@length,varlist));
	pformt_single = ['%-' num2str(maxl+5) '.' num2str(maxl) 's'];
	pformt = ['  ' repmat(pformt_single,1,ncols) '\n'];
	
	n = numel(varlist);
	indResort = reshape(1:ceil(n/ncols)*ncols, [], ncols)';
	
	fprintf(1,pformt,varlist{indResort(indResort<=n)})
	if mod(numel(varlist),ncols)
		fprintf(1,'\n')
	end
end

function printVerbose(varlist)
	% This prints vars sorted by Bytes.
	
	expr = sprintf('evalc(''whos %s'')', sprintf('%s ',varlist{:}));
	z = evalin('base',expr);
	if isempty(z)
		return
	end
	zz = strsplit(z,'\n');
	ne = find(~cellfun('isempty',zz));
	
	header = zz{1};
	vars = zz(ne(2:end));
	indEnd = regexp(header, 'Bytes', 'end');
	
	t = cellfun(@(x)x(1:indEnd), vars, 'UniformOutput', false);
	tt = regexp(t, '\d+$', 'match');
	bytes = str2double([tt{:}]);
	[~,ia] = sort(bytes);
	% Empty lines that 'whos' prints are intentially left of of this output.
	%	This also prints header in orange.
	fprintf(1,'[\b%s]\b\n',header)
	fprintf(1,'%s\n', vars{ia})
	
	% This prints vars by alphabetical order of Name.
% 	expr = sprintf('whos %s', sprintf('%s ',varlist{:}));
% 	evalin('base',expr)
% 	fprintf(1,'\b')
end

function clearThese = getUnprotectedVars(protectedVars)
	%	current_vars = evalin('caller','who')';
	current_vars = evalin('base','who')';
	clearThese = current_vars( ...
		~ismember(current_vars,protectedVars) );
end
