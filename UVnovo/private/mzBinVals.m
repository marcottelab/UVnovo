function [mzrankPreds, colNames] = mzBinVals(specIn, nodesFwdRev, scan2node, mzOffsets, mzThresh, accumArgs)
% specIn <n x 1 struct> array of spectra peaks.
%	fields: '.mz', '.(accumArgs{m,1})'
%	There can be multiple different fields dynamically referred to in accumArgs.
% % specNodes <n x 1 cell> float array of bin centers for each spectrum
% nodesFwdRev <n x 1 cell> float array [b x 2] mass nodes for each spectrum
%	[nterm, cterm] charge 1 masses.
%	Each bin is centered at a specNode + mzOffset
% scan2node <n x 1 int> map from specIn to nodesFwdRev
% mzOffsets <z x 1 struct> of mass offsets from N & C term
%	Index (z) matches precursor charge.
%	fields: '.nterm', '.cterm'
% mzThresh <1 x 2 float> mz threshold, lo & hi, for m/z bins
% accumArgs <m x 3 cell> 2nd, 4th, 5th args to accumarray builtin function
%	{VAL, FUN, FILLVAL}
%	{m, 1} <str> is a field name in specIn(n).peaks.(accumArgs{m,2})
%		Each of these fields must match dimensions of specIn(n).peaks.mz.
%	{m, 2} <function handle> is an accumulation function
%	{m, 3} <same class as FUN return> value to put in empty bins
%	ex. {'rank', @min, nan} return the min ranked peak for each bin
% 
% mzrankPreds <n x m cell> of [npeaks x nOffsets]
%	columns ordered as [nterm_z1, cterm_z1, nterm_z2, cterm_z2, ...]
% 
% @TODO this doc stuff is partly not current any more.
% @TODO rename some of the vars that are specific for predvec stuff
% 
% @TODO Has much slow. Should have less slow.
% This function is MUCH faster when nodesFwdRev contains only the unique set &
%	not a 1:1 mapping to specIn. I'm not sure why there's so much difference.
%	Memory & caching effects?
%	Ex. Speed comparison for 4 parallel workers:
%		numel(nodesFwdRev)==7911 and matches size specIn -- 50% slower
%		numel(nodesFwdRev)==1245 and maps to 7911 specIn
% 
% !!! replace all the indexing.range stuff with histc & fancy indexing & vectorization
%	See mzBinVals_devhistc.m

UNIT_G = CONSTS.unit_g;
nscans = numel(specIn);
naccum = size(accumArgs,1);
npcharge = numel(mzOffsets);

%	If any of specIn, accumArgs, or mzOffsets are empty, return.
if ~nscans || ~naccum || ...
		all(arrayfun(@(i)all(~structfun(@numel,mzOffsets(i))),1:npcharge))
	mzrankPreds = cell(nscans,0);
	colNames = {};
	return
end

% % specFields = fieldnames(specIn(1).peaks);
specFields = fieldnames(specIn);
for i = 1:naccum
	f = accumArgs{i,1};
	assert(ismember(f,specFields), exceptionID('specFieldMissing'), ...
		sprintf('Field ''%s'' not present in struct <ARG1>', f));
end
%%

nBinSets = numel(nodesFwdRev);

binSets = cell(nBinSets,1);
ib = cell(nBinSets,1);
pVecDims = zeros(nBinSets,2);
c_binCenters = cell(npcharge,2);
for i = 1:nBinSets
	nodes = nodesFwdRev{i}(:,1);
	nodesRev = nodesFwdRev{i}(:,2);
	for z = 1:npcharge
		% Proton mass diff when converting mh->m/z. Slight approximation gets 
		%	worse when unit_g diverges farther from mProton.
		e = (z-1)/z;
		c_binCenters{z, 1} = bsxfun(@plus, nodes/z + e, mzOffsets(z).nterm);
		c_binCenters{z, 2} = bsxfun(@plus, nodesRev/z + e, mzOffsets(z).cterm);
	end
	t_binCenters = [c_binCenters{:}];
	[t_uniInt, ~, ib{i}] = unique(t_binCenters);
	t_uniReal = t_uniInt*UNIT_G;
	binSets{i} = [t_uniReal + mzThresh(1), t_uniReal + mzThresh(2)];
	pVecDims(i,:) = size(t_binCenters);
end

binSetsExpand = binSets(scan2node);
ibExpand = ib(scan2node);
pVecDimsExpand = pVecDims(scan2node,:);

mzrankPreds = cell(nscans,naccum);
parfor n = 1:nscans
% 	n_binCenter = scan2node(n);
% 	binSet = binSets{n_binCenter};
% 	accSz = [size(binSet,1),1];
	accSz = [size(binSetsExpand{n},1),1];
	[~, indSpec, indBin] = ranges(specIn(n).mz, binSetsExpand{n});
	for i = 1:naccum
		accVals = specIn(n).(accumArgs{i,1});
		accFun  = accumArgs{i,2};
		accFill = accumArgs{i,3};
		
		a = accumarray(indBin, accVals(indSpec), accSz, accFun);
		if ~isequal(accFill, 0)
			a(a==0) = accFill; % !!!!!! This is faster than providing accumarray a different fillval, but it will cause problems if accumarray legitimately places a 0 !!!!!!
		end
		mzrankPreds{n,i} = reshape(a(ibExpand{n}), pVecDimsExpand(n,:));
	end
end

if nargout > 1
	colNames = convertPredNames_mzOffsets(mzOffsets);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [rangeInds, indXX, indBin] = ranges(xx, boundSets)
	% @TODO finish documenting or get rid of completely.
	% @TODO rename the input & output vars to something more useful.
	% 
	% INPUT
	%	xx <vector> sorted numeric vector
	%	boundSets <array> n x 2 num array. Each numeric pair represents lo (n,1) 
	%		and hi (n,2) range bounds.
	% OUTPUT
	%	rangeInds <cell>
	%	indXX <vector> t x 1 numeric index into xx
	%	indBin <vector> t x 1 numeric index of range bin for elmnts indexed in indXX
	% 
	% indXX, indBin are used to map between elements of xx and ranges in boundSets.
	%	These have length [t] equivalent to number of elements in xx that fall into 
	%		ranges in boundSets, allowing duplicates and overlapping bins.
	% binMapping usage examples
	%	Accumulate elements falling within bin ranges & apply a reducing funct:
	%		yy = accumarray(indBin, xx(indXX)); % sums members in each bin
	%		yy = accumarray(indBin, xx(indXX), size(xx), @min, nan);
	%		yy = accumarray(indBin, zz(indXX), size(xx)); % zz - vec same size as xx

	if isempty(xx)
		rangeInds = zeros(0,2);
		indXX = zeros(0,1);
		indBin = zeros(0,1);
		return
	end

	if isempty(boundSets)
		rangeInds = zeros(0,2);
		indXX = zeros(0,1);
		indBin = zeros(0,1);
		return
	end

	assert(issorted(xx))

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	indbp2all = boundSets(:,1) < boundSets(:,2);
	%	toss out pairs where window width would be <= 0
	bpairs = boundSets(indbp2all,:);

	%	sort all values in bpairs. Keep track of original index.
	[bb, ibp] = sort(bpairs(:));
	%	indices back to bound pair matrix
	islowind = ibp <= size(bpairs,1); % starting b values (vs ending)

	%%% iterate through xx and bb
	nbb = length(bb);
	nxx = length(xx);
	ix = 1;
	ib = 1;
	xb = zeros(nbb,1);

	while ib <= nbb && bb(ib) < xx(1)
		ib = ib+1;
	end

	% The linear search in this while block is better when bins are fairly close.
	% Otherwise the partly-binary search, commented at bottom, is probably faster.
	% A C function would be ideal, or getting rid of this altogether.
	while ib <= nbb
		while bb(ib) > xx(ix) && ix < nxx % this is slow but for loop is much slower
			ix = ix+1;
		end
		if bb(ib) > xx(ix)
			xb(ib:end) = ix + islowind(ib:end);
			break
		end
		if islowind(ib)
			xb(ib) = ix;
		else
			xb(ib) = ix-1;
		end
		ib = ib+1;
	end

	xxbb = zeros(size(bpairs));
	xxbb(ibp) = xb;

	rangeInds = zeros(size(boundSets));
	rangeInds(indbp2all,:) = xxbb;
	rangeInds(rangeInds(:,1)==0,1) = 1;
	rangeInds(rangeInds(:,1)>rangeInds(:,2),:) = 0;

	if nargout > 1
		maxBinMems = max(diff(rangeInds,[],2))+1;
		% Storing counterfill in a persistent isn't any faster.
		counterfill = cell(1,maxBinMems);
		for j = 1:maxBinMems
			counterfill{j} = [0:j-1];
		end

		indNE = find(rangeInds(:,1)>0);
		rb = rangeInds(indNE,:);
		d = diff(rb,[],2);
		t = [counterfill{d+1}]';
		ind_r = cumsum(~t);

		indXX = rb(ind_r,1) + t;
		indBin = indNE(ind_r);
	end

	% %%% Partly binary search start
	% a0 = 4;
	% while ib <= nbb
	% 	bb_ib = bb(ib);
	% 	a=a0;
	% 	while ix < nxx && bb_ib > xx(ix) % binary search, start with interval 'a'
	% 		a = a*2;
	% 		ix = ix+a;
	% 	end
	% 	if ix > nxx && bb_ib > xx(nxx) % then we're done
	% 		ix = nxx;
	% 	elseif a<=a0*4
	% 		if a>a0
	% 			ix=ix-a;
	% 		end
	% 		while bb_ib > xx(ix) && ix < nxx % back to the linear search...
	% 			ix = ix+1;
	% 		end
	% 	else
	% 		a=a/2;
	% 		ix = ix-a;
	% 		while ix > nxx || (bb_ib <= xx(ix) && a>1)
	% 			a=a/2;
	% 			ix = ix-a;
	% 		end
	% 		if a>16 % then narrow linear search to correct 1/16th of current range
	% 			a = a/8;
	% % 			a=8;
	% 			ix = ix+a;
	% 			while ix < nxx && bb_ib > xx(ix) % binary search, start with interval 'a'
	% 				ix = ix+a;
	% 			end
	% 			ix=ix-a;
	% 		end
	% 		while bb_ib > xx(ix) && ix < nxx % back to the linear search...
	% 			ix = ix+1;
	% 		end
	% 	end
	% 	if ix == nxx && bb_ib > xx(ix)
	%  		xb(ib:end) = ix + islowind(ib:end);
	% 		break
	% 	end
	% 	if islowind(ib)
	% 		xb(ib) = ix;
	% 	else
	% 		xb(ib) = ix-1;
	% 	end
	% 	ib = ib+1;
	% end
	% %%% end of binary index searching
end
