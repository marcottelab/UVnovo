function vars = loadSer(filename, asStruct)
% LOADSER(filename)
%	Load serialize data (that was saved with io.saveSer).
%	This is MUCH faster and more space-efficient for non-numeric data!
% 
% 
%	This WILL overwrite vars in calling workspace if of same name.
%	Use option 'asStruct' if that may be a problem.
% 
% See also IO.SAVESER


if ~exist('asStruct','var') || isempty(asStruct), asStruct = false; end

assert(logical(exist(filename,'file')), ...
	'Unable to read ''%s'': no such file or directory.', filename)

builtin('load',filename);

assert(logical(exist('bs','var')), ...
	'File to load must originate from io.saveSer() or match its format.')

varStruct = getArrayFromByteStream(bs);
if asStruct
	vars = varStruct;
	return
end

vars = fieldnames(varStruct);
for i = 1:numel(vars)
	f = vars{i};
	assignin('caller', f, varStruct.(f))
end

clear varStruct
if ~nargout
	clear vars
end