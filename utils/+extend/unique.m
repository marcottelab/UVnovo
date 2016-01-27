
function structOut = unique(arrayIn,uniqueFlag,calcFieldsIn)
% extend.unique() add some functionality I commonly use with 'unique'.
%	Ex. tabulating counts, mapping between unique and original, etc.
% 
% INPUT ARGS
%	arrayIn
% 	uniqueFlag: <str> arg passed to base 'unique'.
%		<'rows', 'first', 'last', 'stable', 'sorted', 'legacy', 'R2012a'>
%	calcFields <cellstr> [optional] fields to calculate & return.
%		<uni, ind_uni2full, ind_full2uni, uni2full, uni2fullsubs, counts>
% OUTPUT
%	structOut
%
% 'stats.tabulate' does some of this. I didn't realize it existed.
%
% @TODO documentation
% @TODO review code. I've probably learned since its implementation.
% @TODO either implement 'calcFieldsIn' option or remove completely.

%%
calcFields={'uni','ind_uni2full','ind_full2uni','uni2full','uni2fullsubs','counts'};
if exist('calcFieldsIn','var') && iscellstr(calcFieldsIn)
	warning('extend_unique_:BrokenCalcFields','This is broken. Fix.')
% 	tf = cellfun(@(x) any(strcmpi( calcFieldsIn,x )), calcFields);
% 	sFields = cell2struct(num2cell(tf), calcFields);
	
	% Calculate all fields until 'calcFieldsIn' option works.
	sFields = cell2struct(num2cell(true(numel(calcFields),1)), calcFields);
else
	sFields = cell2struct(num2cell(true(numel(calcFields),1)), calcFields);
end

if exist('uniqueFlag','var')
	if ischar(uniqueFlag)
		if ismember(lower(uniqueFlag),...
				{'rows','first','last','stable','sorted','legacy','R2012a'})
			[uni, ia, ic] = unique(arrayIn,lower(uniqueFlag));
		else
			warning('extend_unique:UniqueFlag_invalid','uniqueFlag is invalid')
			[uni, ia, ic] = unique(arrayIn);
		end
	else
		warning('extend_unique:UniqueFlag_badType','uniqueFlag must be a string')
	end
else
	uniqueFlag = '';
	[uni, ia, ic] = unique(arrayIn);
end

nUni = numel(uni);

if sFields.uni  % sorted set of unique elements
	structOut.uni = uni;
end
if sFields.ind_uni2full  % index of each unique to first instance in arrayIn
	structOut.ind_uni2full = ia;
end
if sFields.ind_full2uni  % index of each element in arrayIn to index in uni
	try
		structOut.ind_full2uni = reshape(ic,size(arrayIn));
	catch
		structOut.ind_full2uni = ic;
	end
end

if sFields.uni2full || sFields.uni2fullsubs
	if strcmpi(uniqueFlag,'rows')
		uni2full = [];
	else
		fullinds = accumarray( ic, 1:size(ic,1), [nUni,1], @(x){x});
		if ischar(uni)
			uni2full = containers.Map(cellstr(uni(:)), fullinds);
		elseif isempty(uni)
			uni2full = containers.Map();
		elseif isa(uni,'categorical')
			uni2full = containers.Map(getlabels(uni), fullinds);
		else
			uni2full = containers.Map(uni, fullinds);
		end
	end
	
	if sFields.uni2full % mapping of each unique element to all inds in arrayIn
		structOut.uni2full = uni2full;
	end
	if  sFields.uni2fullsubs % mapping of each unique element to subs in arrayIn
		if nnz(size(arrayIn)>1)==1 || isempty(uni2full)
			structOut.uni2fullsubs = uni2full;
		else
			fs = cell(nUni,1);
			for i = 1:nUni
				[fs{i}(:,1), fs{i}(:,2)] = ind2sub(size(arrayIn),fullinds{i});
			end
			uni2fullsubs = containers.Map(uni, fs);
			structOut.uni2fullsubs = uni2fullsubs;
		end
	end
end
if sFields.counts  % counts of each unique element
	structOut.counts = accumarray(ic,1);
end

