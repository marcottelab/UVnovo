function [mzOffsets, offsetNames, nonoffsetNames] = convertPredNames_mzOffsets(predNamesORmzOffsets)
% convertPredNames_mzOffsets: convert between mzOffsets struct & predictor names
%	This assumes predNames are ordered correctly.
% mzOffsets <z x 1 struct>
% offsetNames <cellstr> predictor names corresponding to values in mzOffsets
% nonoffsetNames <cellstr> predictor names that were not parseable to mzOffsets
% 
% This doesn't yet handle case where predNames indicates multiple spectra types.
% 
%%%%%%%%%%%%%%%
% @TODO work on documentation and making this function more idiomatic of UVnovo.


if iscellstr(predNamesORmzOffsets)
	convertTo = 'mzOffsets';
	predNames = predNamesORmzOffsets;
elseif isstruct(predNamesORmzOffsets) && ...
		any(isfield(predNamesORmzOffsets,{'nterm','cterm'}))
	convertTo = 'predNames';
	mzOffsets = predNamesORmzOffsets;
else
	error(exceptionID('unknownInputType'), ...
		'Input must be a cellstr of predictor names or mzOffsets struct')
end

%%
% varargout = cell(nargout);
switch convertTo
	case 'predNames'
		% create predictor names from mzOffsets struct
		npcharge = numel(mzOffsets);
		t_colNames = cell(npcharge,2);
		for z = 1:npcharge
			if isfield(mzOffsets,'nterm')
				expr_nz = sprintf('nTerm_z%d_%%gmz', z);
				t_colNames{z, 1} = sprintfc(expr_nz, mzOffsets(z).nterm);
			end
			if isfield(mzOffsets,'cterm')
				expr_cz = sprintf('cTerm_z%d_%%gmz', z);
				t_colNames{z, 2} = sprintfc(expr_cz, mzOffsets(z).cterm);
			end
		end
		offsetNames = [t_colNames{:}];
		nonoffsetNames = cell(1,0);
		
	case 'mzOffsets'
		% create mzOffsets struct from predictor names
		mzOffsets = struct('nterm',[], 'cterm',[]);
		
		expr = '([nc])Term_z(\d+)_([-\+]?\d*\.?\d*)mz(.*)';
		c_toks = regexp(predNames, expr, 'tokens');
		indHits = ~cellfun('isempty', c_toks);
		if ~any(indHits) % then return
			warning(exceptionID('noValidOffsetPredNames'), ...
				'No predictor names recognized as offset predictors.')
			nonoffsetNames = predNames;
			offsetNames = cell(1,0);
			return
		end
		
		t = cat(1,c_toks{indHits});
		tokens = cat(1,t{:});
		
		ncTerm = [tokens{:,1}]';
		pcharges = str2double(tokens(:,2));
		offsets = str2double(tokens(:,3));
		
		offsetNamesIn = predNames(indHits);
		specset = tokens(:,4); % ex. this could represent activation ('UVPD')
		if any(~cellfun('isempty',specset))
			warning(exceptionID('TODO'),'IMPLEMENT SPECSET DIFFERENTIATION')
			
			offsetNames_n = regexprep(offsetNamesIn, ...
				'([-\+]?\d*\.?\d*mz).*', '$1');
		else
			offsetNames_n = offsetNamesIn;
		end
		
		z2pred = accumarray(pcharges,1:numel(pcharges),[],@(x){x});
		for z = numel(z2pred):-1:1
			iz = sort(z2pred{z});
			if isempty(iz), continue, end
			inz = ncTerm(iz) == 'n';
			if ~any(inz)
				mzOffsets(z).nterm = zeros(1,0);
			else
				mzOffsets(z).nterm = offsets(iz(inz))';
			end
			if all(inz)
				mzOffsets(z).cterm = zeros(1,0);
			else
				mzOffsets(z).cterm = offsets(iz(~inz))';
			end
		end
		
		[~,offsetNamesNew,~] = convertPredNames_mzOffsets(mzOffsets);
		if ~isequal(offsetNames_n(:),offsetNamesNew(:))
			[~,ia] = ismember(offsetNamesNew, offsetNames_n);
			if ~issorted(ia)
				warning(exceptionID('inputPredVarsNotSorted'), ...
					'Input predictor order does not match output.')
			else
				error(exceptionID('predNameMismatch'), ...
					'Inconsistent predictor variables! This should never occur.')
			end
			offsetNames = offsetNamesIn(ia);
		else
			offsetNames = offsetNamesIn;
		end
		nonoffsetNames = predNames(~indHits);
end

