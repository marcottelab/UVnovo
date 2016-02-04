function [h, importanceTable] = plotPredictorImportance(varargin)
% plotPredictorImportance(ensTB, name-val pairs)
%	ensTB: TreeBagger ensemble, built with option ('oobvarimp', 'on')
% plotPredictorImportance(H, ensTB, name-val pairs)
%	H: figure or axes handle
% name-val pairs:
%	'score' - score metric to plot <'OOBPermutedVarDeltaError' [default], 
%		'DeltaCritDecisionSplit', 'OOBPermutedVarDeltaMeanMargin',
%		'OOBPermutedVarCountRaiseMargin'>
%		The last 2 scores metrics are nearly identical to the default,
%			from what I've seen.
% h: handle
% importanceTable: <cell (n x 3)> predictor name, rank, importance score
% 
% This func could use better documentation.


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Parse arguments & initialize figure.

if nargin > 1 && ~ischar(varargin{2})
	try % first arg can be a figure or axes handle, or it can be a scalar int.
		h = varargin{1};
		if ishghandle(h) && strcmpi(get(h,'Type'),'axes')
			axes(h); % set current axes to provided axes handle
			cla
		else
			figure(h);
			clf
		end
		hold on
		cla
		ensTB = varargin{2};
		assert(isa(ensTB, 'TreeBagger') || isa(ensTB, 'CompactTreeBagger'), ...
			'MS:RF:plotPredictorImportance:badInput', ...
			'second argument must be of class TreeBagger when passing figure handle as first argument')
	catch me
		rethrow(me)
	end
	varargin(1:2) = [];
else
	h = figure; hold on
	ensTB = varargin{1};
	assert(isa(ensTB, 'TreeBagger') || isa(ensTB, 'CompactTreeBagger'), ...
		'MS:RF:plotPredictorImportance:badInput', ...
		'first argument must be a TreeBagger object or figure handle')
	varargin(1) = [];
end


% Optional parameters. Defaults:
scoreVar = 'OOBPermutedVarDeltaError';
transformYScores = false;
y_label = 'Importance';

% Run through remaining varargin for name-val user params.
while numel(varargin)
	optname = varargin{1};
	optarg = varargin{2};
	switch optname
		case 'score'
			if ismember(optarg, {
					'OOBPermutedVarDeltaError'; % default
					'DeltaCritDecisionSplit';
					'OOBPermutedVarDeltaMeanMargin';
					'OOBPermutedVarCountRaiseMargin'} )
				scoreVar = optarg;
			end
		case 'yScaleFactor'
			assert(isscalar(optarg),[mfilename, ':BadYScaleFactorArg'], ...
				'yScaleFactor must be a scalar')
			yScaleFactor = optarg;
			if yScaleFactor~=1
				yTransform = @(y)y.^yScaleFactor;
				transformYScores = true;
				y_label = sprintf('Importance scaled by power(%g)', yScaleFactor);
			end
		case 'yTransform'
			assert(isa(optarg,'function_handle'), [mfilename, ':BadYTransformArg'], ...
				'yTransform arg must be a function handle')
			yTransform = optarg;
			transformYScores = true;
			y_label = 'Importance (transformed score)';
		otherwise
			warning([mfilename, ':UnknownParam'], ...
				'Unknown parameter name ''%s''', optname)
	end
	varargin(1:2) = [];
end

xlabel('Predictor')
ylabel(y_label)
title(sprintf('Predictor Importance, using %s', scoreVar))


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Get predictor importance.

rawScores = ensTB.(scoreVar); % this is slow

if transformYScores
	try
		predScores = yTransform(rawScores);
	catch
		error([mfilename, ':BadYTransformFunction'], ...
			'yTransform function must act on a vector of doubles')
	end
else
	predScores = rawScores;
end


[~, ia] = sort(rawScores,'descend');
predScoresSorted = predScores(ia);
nPreds = numel(ia);

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Plot predictor importance.

hBar = bar(predScores,'k');
hPlot = plot(predScoresSorted(end:-1:1), 'color',[1, 0.41015625, 0.703125], ...
	'LineSmoothing','on', 'LineWidth',2);

% Assign user data for custom datatips.
prank(ia,1) = 1:nPreds;

if transformYScores
	dataFormat = {'%s','rank: %4d','raw score: %6g','scaled: %6f'};
	barData = [ensTB.VarNames', num2cell(prank), ...
		num2cell(rawScores)', num2cell(predScores)'];
	plotData = flipud([ensTB.VarNames(ia)', num2cell(1:nPreds)', ...
		num2cell(rawScores(ia))', num2cell(predScoresSorted)']);
else
	dataFormat = {'%s','rank: %4d','score: %6g'};
	barData = [ensTB.VarNames', num2cell(prank), num2cell(predScores)'];
	plotData = flipud([ensTB.VarNames(ia)', num2cell(1:nPreds)', ...
		num2cell(predScoresSorted)']);
end
set(hBar, 'UserData', struct('data', {barData},  'format', {dataFormat}));
set(hPlot,'UserData', struct('data', {plotData}, 'format', {dataFormat}));

% Set custom data cursor callback.
dcm_obj = datacursormode(gcf);
set(dcm_obj,'UpdateFcn',@UserData_cursor, 'Enable','on')

% Make figure a bit nicer.
set(gcf,'Renderer','OpenGL')
xlim([.5, nPreds+.5])

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Finish

% This looks better than 'vargout' and accomplishes same thing.
if nargout == 0
	clear h
elseif nargout >1
	importanceTable = barData(:,1:3);
end
