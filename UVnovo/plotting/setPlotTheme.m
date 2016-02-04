function setPlotTheme(h, themeNew, themeCurr)
% SETPLOTTHEME change the thematic style of <h> and all its children.
%	This currently only changes plot colors. It can be extended to other
%	properties, too.
%	The implementation is pretty basic. It could be made more flexible, robust,
%	etc with a little work but suits my needs as it is.
% 
% INPUT ARGS:
%	h: graphics handle (fig or axes). This is typically the root figure handle.
%	themeNew: desired theme
%	themeCurr: current theme. This can be detected from h.UserData if function
%		was already executed for the handle object.
% 
%	Implemented themes:
%		DEFAULT: standard, ugly Matlab plots.
%		DARK: dark. yup. Not good for printing.
% 

%% Themes
%	matlab default
theme.default.text = [0, 0, 0];
theme.default.background = [.8, .8, .8];
theme.default.axes = [1, 1, 1];
%	dark
theme.dark.text = [.7, .7, .7];
theme.dark.background = [.168, .168, .168];
theme.dark.axes = [.108, .108, .108];

themes = fieldnames(theme);

%% Init input arg vars

if ~exist('h','var') || ~isscalar(h) || ~ishghandle(h)
	warning('setPlotTheme:nonscalarArg1','Arg 1 must be scalar graphics handle')
	return
end
if ~exist('themeNew','var') || isempty(themeNew), themeNew = 'dark'; end

% Properties can be accessed like a struct when h is explicitly a handle.
hh = handle(h); % hh and h reference the same data
if ~exist('themeCurr','var') || isempty(themeCurr)
	% isfield works even on non-structs (returning false)
	if isfield(hh.UserData, 'plotTheme') ...
			&& ischar(hh.UserData.plotTheme) ...
			&& any(strcmpi(hh.UserData.plotTheme, themes))
		themeCurr = hh.UserData.plotTheme;
	else
		themeCurr = 'default';
	end
end


%%% Main %%%

colorPropNames = {'Color','XColor','YColor','ZColor','EdgeColor','TextColor'};

fieldsCurr = fieldnames(theme.(themeCurr));
fieldsNew =  fieldnames(theme.(themeNew));
[a, ia] = ismember(fieldsCurr, fieldsNew);
aa = struct2cell(theme.(themeCurr));
bb = struct2cell(theme.(themeNew));
colorCurrNew = [aa(a), bb(ia)];

nc = size(colorCurrNew, 1);
hChange = cell(nc, numel(colorPropNames));

for p = numel(colorPropNames):-1:1
	prop = colorPropNames{p};
	for c = nc:-1:1
		hChange{c,p} = findall(h, prop, colorCurrNew{c,1});
	end
end

for p = 1:numel(colorPropNames)
	prop = colorPropNames{p};
	for c = 1:nc
		hColor = hChange{c,p};
		set(hColor, prop, colorCurrNew{c,2})
	end
end

