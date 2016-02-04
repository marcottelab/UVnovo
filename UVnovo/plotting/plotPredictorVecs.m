function h = plotPredictorVecs(pv, varNames, varargin)
% PLOTPREDICTORVECS plots ensemble predictor variables with custom data tips.
%	Useful for data introspection.
% 
%	PV <array [nObs x nVars]> numeric matrix of predictor vars.
%	VARNAMES <cellstr {1 x nVars}> predictor names.
%	VARARGIN <see initParams> optional params, as described below for paramsDef.
% 
% This could be fleshed out a little.

paramsDef = struct(  ...
	'clim', [-1 1],  ... % value range for colorbar scaling.
	'ignoreVal',nan, ... % This val is ignored when summing across observations.
	'varInds', [0],  ... % Only plot subset of pv vars. <0 [default]> plot all.
	'subsample', 1,  ... % Frac [x <= 1] or max number [x > 1] of obs to plot.
	'colorMids', 4   ... % controls colormap gradient. Higher -> extended lows.
	);
% if ~exist('paramsIn','var'), paramsIn = []; end
params = initParams(paramsDef, varargin{:});

if ~exist('varNames','var') || isempty(varNames)
	varNames = sprintfc('var %g', 1:size(pv,2));
end

% Randomly subsample observations.
if params.subsample <= 1
	pvSubInds = find( rand(size(pv,1),1) <= params.subsample );
else
	nObs = min(size(pv,1), ceil(params.subsample));
	pvSubInds = sort( randperm(size(pv,1), nObs) );
end

if isnumeric(params.varInds) && ~isequal(0, params.varInds)
	% Retain specified vars.
	pv = pv(pvSubInds, params.varInds);
	varNames = varNames(params.varInds);
else
	pv = pv(pvSubInds,:);
end

[nObs, nVars] = size(pv);

h = figure;

% Plot 2D image of predictors.
subplot(4,1,1:3)
hAx1 = imagesc(pv, params.clim);
set(hAx1, 'UserData', struct( ... % data and formatting for custom cursor
	'data', {{
				repmat(varNames, nObs, 1);
				pv;
				repmat(pvSubInds, 1, nVars);
			}}', ...
	'format', {{'%s','value: %g','obs: %g'}}));
dcm_obj = datacursormode(gcf);
set(dcm_obj,'UpdateFcn',@UserData_cursor, 'Enable','on')

% Plot sum along y dimension.
subplot(4,1,4)
pv_tosum = pv;
pv_tosum(pv==params.ignoreVal) = 0;
hAx2 = plot(sum(pv_tosum));
set(hAx2, 'UserData',struct('data',{{varNames}}, 'format',{{'%s'}}));
dcm_obj = datacursormode(gcf);
set(dcm_obj,'UpdateFcn',@UserData_cursor, 'Enable','on')

% Set figure properties.
cm = squeeze(real2rgb(1:64, ...
	[[0 0 0; .4 .33 .33; .78 .78 .72; 1 1 1], (4:-1:1)'.^params.colorMids], []));
colormap(cm)

setPlotTheme(h, 'dark')

linkaxes(cell2mat(get([hAx1,hAx2],'Parent')),'x');
zoom xon

if nargout==0
	clear h
end
