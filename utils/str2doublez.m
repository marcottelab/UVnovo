function arrayOut = str2doublez(cstrs)
% STR2DOUBLEZ Convert strings in cellstr arary CSTRS to double precision value.
%   Each string in CSTRS should be an ASCII character representation of a real
%   value. It may only contain digits, a decimal point, and a leading + or -. It
%   may additionally be 'NaN' to return NaN.
% 
% This is not as flexible as MATLAB's builtin str2double, but can be >40x faster
%	on large arrays.
%
%   See also STR2NUM, STR2DOUBLE.

if ~iscellstr(cstrs)
	arrayOut = str2double(cstrs);
else
	%	Pad strings so that each num in char array is delimited for sscanf.
	[~,ia] = max(cellfun('length',cstrs));
	cstrs{ia}(end+1) = ' ';
	arrayOut = reshape(sscanf( char(cstrs)', '%f'), size(cstrs));
end

