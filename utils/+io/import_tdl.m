
function [data, header, headerVarNames, columnInds] = import_tdl(filenameIn, hasHeader)
% import_tdl: load a tab-delimited file that has uniform # columns
%	headerVarNamed uses custom translation of header into valid variable names
%	columnInds is index struct into column fields
%	pretty basic, could be expanded
%		-check if header is present
% 
% @TODO deprecate the need to use this and use a builtin function instead.

if ~exist('hasHeader','var') || isempty(hasHeader), hasHeader = true; end

% params.uniqueFlag = 'stable'; % flag (string) for base unique function

fidIn = fopen(filenameIn);
z_cleanupObj = onCleanup(@() fclose('all')); % just in case

fcontent=fread(fidIn,'*char')';
fclose(fidIn);

fcontent = regexprep(fcontent,'\r?\n?$','');
t_split = regexp(regexp(fcontent,'\r?\n?','split'),'\t','split');
ts_sizes = cellfun('size',t_split,2);
if any(ts_sizes>1)
	t_split(ts_sizes==1) = [];
end


if hasHeader
	header = t_split{1};
	data = cat(1,t_split{2:end});
else
	header = {};
	data = cat(1,t_split{1:end});
end

if nargout>2
	expRep = { '[: ; , \/ \+ -]','_'; '\.',''; '#','n'; '__','_' };
	t_header = header;
	for nrep = 1:size(expRep,1)
		t_header = regexprep(t_header,expRep{nrep,1},expRep{nrep,2});
	end
	headerVarNames = genvarname(t_header);
	if nargout>3
		columnInds = io.cinds(headerVarNames,'_');
	end
end
