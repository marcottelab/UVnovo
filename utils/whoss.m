function whoss(varargin)
% WHOSS prints output of whos sorted by increasing bytes.
%	Empty lines that 'whos' prints are intentially left of of this output.
%	Also, the header prints in orange.
% 
%	See also WHOS


% Parse the printed output of WHOS rather than its return struct so that the
%	'Size' and 'Attributes' columns are retained.
expr = sprintf('evalc(''whos %s'')', sprintf('%s ',varargin{:}));
z = evalin('base',expr);
if isempty(z)
	return
end
zz = strsplit(z,'\n');
ine = find(~cellfun('isempty',zz));

header = zz{1};
vars = zz(ine(2:end));
indEnd = regexp(header, 'Bytes', 'end');

t = cellfun(@(x)x(1:indEnd), vars, 'UniformOutput', false);
tt = regexp(t, '\d+$', 'match');
bytes = str2double([tt{:}]);
[~,ia] = sort(bytes);

fprintf(1,'[\b%s]\b\n',header)
fprintf(1,'%s\n', vars{ia})




