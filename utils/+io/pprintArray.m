function formattedOutput = pprintArray(arrayIn, colNames, rowNames, paramsIn)
% Pretty Print array
% io.pprintArray
% io.pprintArray(arrayIn, colNames, rowNames, struct('format','%d'))
% io.pprintArray(arrayIn, colNames, rowNames, {'format','%d', 'printOrange',0})
% 
% Only prints cells containing single values or empty cells.
%	Modify in future to support nested cells/arrays.
% 
% ARGS
%	colNames: <cell, char, logical, []>
% 
% no real error checking yet
%
% @TODO document
% @TODO structure things better (like the colNames stuff)

%% Parameters
paramsDef = struct( ...
	'indent', '', ...
	'fillna', '', ...
	'format', '%.9g', ...
	'delim', '', ... if provided, place single delim btwn each field
	... @TODO 'wrapCol', 0, ... wrap table at specified column (only if no delim)
	'strtrim', '', ... trim whitespace <default: [true if <delim> else false]>
	'printOrange', true ... print col & row names in orange
	);

if ~exist('paramsIn','var'), paramsIn = []; end
params = initParams(paramsDef, paramsIn);

%%%
if ~exist('colNames','var') || isempty(colNames) || ...
		( islogical(colNames) && colNames == false )
	printColNames = false;
	colNames = '';
else
	printColNames = true;
end

if ~exist('rowNames','var') || isempty(rowNames)
	printRowNames = false;
else
	printRowNames = true;
end

%%% resolve mutually exclusive params
if isempty(params.strtrim) && ~isempty(params.delim)
	params.strtrim = true;
end

if nargout
	%	Then don't print table to command window.
	printArray = false;
	printOrange = false;
else
	printArray = true;
	printOrange = params.printOrange;
end

indent = params.indent;

%% Main

%	Attempt to convert to cellstr if not already.
if isnumeric(arrayIn)
	arrayIn = num2cell(arrayIn);
end

[arrayHeight, arrayWidth] = size(arrayIn);

if printRowNames
	cell_00 = {''};
	rowNames = rowNames(:);
	
	% special cases: user can provide character or numeric arrays to label rows
	if ischar(rowNames)
		rowNames = cellstr(rowNames);
	elseif isnumeric(rowNames)
		rowNames = num2cell(rowNames);
	else
		assert(iscellstr(rowNames))
	end
	
	switch numel(rowNames)
		case 1
			rowNames = repmat(rowNames, arrayHeight, 1);
		case arrayHeight
			
		case arrayHeight + 1
			cell_00 = rowNames(1);
			rowNames = rowNames(2:end);
		otherwise
			warning('io:ppringArray:badvar_rowNames',...
				'var ''rowNames'' doesn''t match arrayIn dimensions')
			t = cell(arrayHeight,1);
			nLabelRows = min( numel(rowNames), arrayHeight);
			t(1:nLabelRows) = rowNames(1:nLabelRows);
			rowNames = t;
	end
else
	rowNames = cell(arrayHeight,0);
	cell_00 =  cell(1,0);
end

% @TODO make code consistent between the row & column name parsing
if printColNames
	
	if isnumeric(colNames)
		colNames = num2cell(colNames);
	elseif ~iscellstr(colNames) % not consistent with row name parsing!
		if ischar(colNames)
			prep = colNames;
		else
			prep = 'col';
		end
		if arrayWidth > 1
			colNames = genvarname(repmat({prep}, 1, arrayWidth), prep);
		else
			colNames = {prep};
		end
	end
	
	colNames = colNames(:)';
	switch numel(colNames)
		case 1
			colNames = repmat(colNames, 1, arrayWidth);
		case arrayWidth
			
		case arrayWidth + 1
			if printRowNames % col takes precedence over row for top left field
				cell_00 =  colNames(1);
			end
			colNames = colNames(2:end);
		otherwise
			warning('io:ppringArray:badvar_colNames',...
				'var ''colNames'' doesn''t match arrayIn dimensions')
			t = cell(arrayWidth,1);
			nLabelCols = min( numel(rowNames), arrayWidth);
			t(1:nLabelCols) = rowNames(1:nLabelCols);
			colNames = t;
	end
else
	colNames = cell(0,arrayWidth);
	if printRowNames
		cell_00 =  cell(0,1);
	else
		cell_00 =  cell(0,0);
	end
end

%%% Compose array to print.
fullArray = [[cell_00, colNames]; [rowNames, arrayIn]];

%	Fill any empty cells.
indEmpties = cellfun('isempty',fullArray);
fullArray(indEmpties) = {params.fillna};

%	Convert all table contents to strings.
indLogicals = cellfun(@islogical,fullArray);
fullArray(indLogicals) = sprintfc('%d',[fullArray{indLogicals}]);
indNums = cellfun(@isnumeric,fullArray);
fullArray(indNums) = sprintfc(params.format,[fullArray{indNums}]);

if params.strtrim
	fullArray = strtrim(fullArray);
end


% @TODO try to get everything below (and everywhere) cleaned up.

%%%	Compose print format strings.
ncols = size(fullArray,2);
colFrmt = cell(1,ncols);

if isempty(params.delim)
	%	Align columns with spaces for pretty printing.
	for c = ncols:-1:1
		colWidthStrs{c} = num2str(max(cellfun('length',fullArray(:,c))));
	end
	colSep = '  '; % padding btwn columns
	
else %	Delimit fields using provided delim char or string.
	colWidthStrs = repmat({''}, 1, ncols);
	assert(ischar(params.delim))
	colSep = params.delim;
end

for c = 1:ncols
	colWidthStr = colWidthStrs{c};
	if c == 1 && printRowNames
		if printOrange
			colFrmt{c} = ['[\b%' colWidthStr 's]\b'];
		else
			colFrmt{c} = ['%' colWidthStr 's'];
		end
	else
		colFrmt{c} = ['%-' colWidthStr 's'];
	end
end
frmtSpec = [strjoin(colFrmt,colSep), '\n'];

if printArray
	strArray_toprint = fullArray';
	if printColNames
		if printOrange
			if printRowNames
				t = regexprep(frmtSpec, {' %', 's '}, {' [\\b%', 's]\\b '});
			else
				t = regexprep(frmtSpec, {'%', 's '}, {'[\\b%', 's]\\b '});
			end
			frmtSpecHeader = [t(1:end-2), ']\b\n'];
		else
			frmtSpecHeader = frmtSpec;
		end
		% could run into problems if 'indent' contains special print chars
		printStr = [sprintf([indent,frmtSpecHeader],strArray_toprint{:,1}), ...
			sprintf([indent,frmtSpec],strArray_toprint{:,2:end})];
	else
		printStr = sprintf([indent,frmtSpec],strArray_toprint{:});
	end
	fprintf(1, printStr)
	
else % Return array as a string.
	t = fullArray';
	formattedOutput = sprintf(frmtSpec,t{:});
end

