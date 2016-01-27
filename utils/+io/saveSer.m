function saveSer(filename, varargin)
% SAVESER(filename, varargin)
%	Serialize data before saving it.
%	This is MUCH faster and more space-efficient for non-numeric data!
% 
% FILENAME: <str> this is passed directly to the builtin SAVE
% VARARGIN: <cellstr> names of vars to save
% 
%	Saved vars are reloaded with io.loadSer
% 
%	ref http://undocumentedmatlab.com/blog/serializing-deserializing-matlab-data
% 
%	!! Anonymous function handles can take FOREVER to save !!
%	This is probably when their state depends on huge vars and has to be saved.
%	Use 'func2str' before saving to remove that overhead. Don't do that if the
%	func is used as a closure and needs to remember state.
% 
%	See also SAVE, IO.LOADSER, IO.CPSER


assert(iscellstr(varargin), 'Argument must contain a string.')

%	Put vars into a single struct object.
nvars = numel(varargin);
for i = 1:nvars
	var = varargin{i};
	try
		v = evalin('caller', var);
	catch
		warning('saveSer:nonexistentVar', '"%s" not present in workspace', var);
		continue
	end
	varStruct.(var) = v;
end

%	Serialize 
bs = getByteStreamFromArray(varStruct); %#ok<NASGU>


%%% Save
%	Performance seems best with v7.3 from a speed/size tradeoff.
%	-v6 followed by 7z is a fair bit smaller but slower (default 7z settings)
%	-v6 followed by zip is marginally smaller (default zip settings)
save(filename, 'bs', '-v7.3');

