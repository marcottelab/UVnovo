function [errorID, status] = exceptionID(mnemonic, stack_location)
% EXCEPTIONID create an MException identifier string, following format of
%	UVnovo:component:mnemonic
% where
%	component: name of function where error is generated
%	mnemonic: message relating to the particular error cause or location
% ARGUMENTS (optional)
%	mnemonic
%	stack_location: <int> position in stack relative to calling function for
%		which error ID should be generated.
% Example: from calling function 'foo'...
%	function foo()
% 	eid = exceptionID('yo',0); % eid == 'UVnovo:foo:yo'

if ~exist('stack_location', 'var') || isempty(stack_location)
	stack_location = 0;
end

s = dbstack(stack_location+1, '-completenames');
if ~isempty(s)
	callingFunc = s.file;
	component = regexprep(callingFunc, ...
		{'\.m$', '^(.*?\\(\+|\@)|.*\\)', '\\(\+|\@)?'}, {'', '', ':'});
	pre = strjoin({'UVnovo',component},':');
else
	pre = 'UVnovo';
end

if exist('mnemonic', 'var') && ~isempty(mnemonic) && ischar(mnemonic)
	errorID = strjoin({pre,mnemonic}, ':');
	status = 1;
else
	errorID = pre;
	status = 0;
end


