function msData = import_spectra(fn_ms2, paramsIn)
% IMPORT_SPECTRA Import spectra from a *.ms2 file, returning the msData
%	structure. This does a little processing, filtering, and standardization of
%	the ms2 data before returning.
% 
% 
% See also LOADMSN.

%%%%%%%%%%%%%%
%%% Parameters

%	This name specifies the primary pmass field.
paramsDef.pmass_field = 'pmass';

%	MS2 peak filtering.
paramsDef.filterPeaks = struct( ...
	'minInten', 5, ...			% remove peaks below this intensity
	'removeZeroInten', true ... % remove peaks with zero intensity
	);

if ~exist('paramsIn','var'), paramsIn = []; end
params = init_params(paramsDef, paramsIn);

%%%%%%%%%%%%%%
%%% Main

fprintf(1, 'Loading file: %s\n', fn_ms2)
t_msData = loadmsn(fn_ms2);

%	Standardize msData field names. @TODO get rid of this.
t_msData.scans = fieldnameStandardize(t_msData.scans);

%	Convert msData fields to formats specified in table:
toconvert = {
	%Field        sscanf format str
	'pmass.*',    '%f';
	'rt',         '%f';
	};
fieldz = fieldnames(t_msData.scans);
for i = 1:size(toconvert,1)
	[fexpr, ftype] = toconvert{i,:};
	convertthese = find( ~cellfun('isempty', regexp(fieldz, fexpr)));
	for j = 1:numel(convertthese)
		convertfield = fieldz{convertthese(j)};
		t = {t_msData.scans.(convertfield)};
		if iscellstr(t)
			t = num2cell(sscanf( strjoin(t), ftype));
			[t_msData.scans.(convertfield)] = t{:};
		end
	end
end

%	Change the default pmass field to alternative. Rename original 'pmass_ORIG'.
if ~isempty(params.pmass_field) && ~strcmp(params.pmass_field, 'pmass')
	[t_msData.scans.pmass_ORIG] = t_msData.scans.pmass;
	[t_msData.scans.pmass] = t_msData.scans.(params.pmass_field);
end

%	Set a unique ID for each scan. Using scan number for now...
[t_msData.scans.uid] = t_msData.scans.nscan;

%	Filter out insignificant peaks.
filteredScans = MSfilterScans({t_msData.scans.spec}, params.filterPeaks);

%	Restructure how spectra peaks are stored.
t = struct;
for i = 1:numel(t_msData.scans)
	mzinten = filteredScans{i};
	t.mz = mzinten(:, 1);
	t.inten = mzinten(:,2);
	%	Rank peaks by decreasing intensity.
	[~,ia] = sort(mzinten(:,2),'descend');
	t.rank = zeros(size(ia));
	t.rank(ia,1) = 1:numel(ia);
	t_msData.scans(i).peaks = t;
	t_msData.scans(i).npeaks = size(mzinten,1);
end
t_msData.scans = rmfield(t_msData.scans,'spec');
t_msData.meta.params = params;


msData = t_msData;

