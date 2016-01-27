function [pathOut, matchingPaths] = genUniquePath(pathIn, suffixDelim)
% genUniquePath() Check if path exists & create new path if so.
%	This path will end with the first positive Int that yields a unique path.
% 
% INPUT ARGS
%	pathIn: <str> directory or file path. May be full path or relative to cd.
%	suffixDelim: <str> delimiter char(s) to use for creating 
% OUTPUT
%	pathOut: <str> unique path
%	matchingPaths: <cellstr> other paths/files matching filename pattern
%%

if ~exist('suffixDelim','var'), suffixDelim = '_'; end

if ~exist(pathIn,'file')
	%	file name is already unique
	pathOut = pathIn;
	matchingPaths = {};
	
else
	[pn, fn, ext] = fileparts(pathIn);
	
	existingPaths = dir(fullfile(pn, [fn, '*', ext]));
	
	%	create reg expression to capture non-matching part of similar file names
	expr_fn =    regexptranslate('escape', fn);
	expr_delim = regexptranslate('escape', suffixDelim);
	expr_ext =   regexptranslate('escape', ext);
	expr = sprintf('^%s%s(.*?)%s$', expr_fn, expr_delim, expr_ext);
	
	%	find matching file names & capture non-matching part
	[t_toks,t_matches] = regexp({existingPaths.name}', expr, 'tokens', 'match');
	indMatches = ~cellfun('isempty',t_matches);
	toks = t_toks(indMatches);
	matches = t_matches(indMatches);
	
	forbiddenSuffs = unique(cellfun(@(x)x{1}{1}, toks, 'UniformOutput', false));
	
	%	get the first unique suffix
	i = 1;
	while ismember(num2str(i), forbiddenSuffs)
		i = i+1;
	end
	
	%	make full file name & path
	uniqueName = sprintf('%s%s%g%s', fn, suffixDelim, i, ext);
	pathOut = fullfile(pn, uniqueName);
	if nargout>1
		matchingPaths = fullfile(pn, unique([{[fn,ext]}, [matches{:}]])');
	end
end

