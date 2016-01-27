

function toCopy = cp(toCopy)
% toCopy = IO.CP(toCopy) Copies output of 'toCopy' to clipboard
%	- useful for surrounding a function or var on command line
%	var passes through if output specified

%%

try
	clipboard('Copy',toCopy);
catch
	clipboard('Copy',...
		strjoin( ...
			cat( 2, ...
				strtrim(regexp(deblank(evalc('disp(toCopy)')),'\n','split')),...
				{''} ),...
			'\n'));
end

if ~nargout
	clear toCopy
end
