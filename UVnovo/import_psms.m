function psmData = import_psms(fn_psms, varargin)
% IMPORT_PSMS import psms from percolator *.xlsx results using load_xlsx_psms.
%	This does additional processing, filtering, and standardization of the data
%	before returning. The excel file is in the Proteome Discover export format.
%	
%	FN_PSMS is the path of a single psm file.
%	VARARGIN optional parameters for load_xlsx_psms.
% 
% See also LOAD_XLSX_PSMS.


%%% Parameters
%	These shouldn't need to change or be visible to end-users.

psmParams = struct( ...
	'rank_limit', 1 ... % maximum rank to include
	);

% This table defines psm columns to import, what to name them in the return
%	struct, and which of potentially many values for a spectrum to retain.
%	Refer to load_xlsx_psms for a more complete description.
%	We set all SELECTION values to 1 to retain only the top PSM of each scan.
psmParams.pFields = {
	% HEADER               NAME        SELECTION
	{'ScanReindexed','FirstScan'} 'nscan' 1; ...
	'Sequence'            'sequence'   1; ...
	'Modifications'       'PTMs'       1; ...
	'XCorr'               'XCorr'      1; ...
	'Rank'                'rank'       1; ...
	'MH0x2B0x5BDa0x5D'    'pmassPD'    1; ... Proteome Discoverer mass estimate
	'Charge'              'pcharge'    1; ...
	};

% @TODO obviate need for 'ScanReindexed' field in the psms file. This is
%	currently used to identify scans after concatenating those from multiple
%	injections. It replaces the run-specific scan number, which may not be
%	unique when combining multiple runs.

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%	Load the psm data.
[t_psmData, psmParams.pFields] = load_xlsx_psms(fn_psms, psmParams, varargin{:});
t_psmData = t_psmData.sample;

%	Annotate sequences with PTM masses.
seqsAnno = annotateSeqs({t_psmData.sequence}, {t_psmData.PTMs});
[t_psmData.seqAnno] = seqsAnno{:};

%	Calculate exact mass of the PSM sequences.
for i = 1:numel(t_psmData)
	t_psmData(i).pmass = MSpepmass(t_psmData(i).seqAnno, [], [], false);
end

%	Compute a hash for each sequence, ignoring PTM variants of a single peptide.
[seqsAA, ~, ic] = unique( upper({t_psmData.sequence}) );
uniseqCounts = accumarray(ic, 1);
nseqs = numel(seqsAA);

h = zeros(nseqs,1);
for i = 1:nseqs
	h(i) = java.lang.String([seqsAA{i}]).hashCode;
end
t = num2cell(h);
[t_psmData.seqHash] = t{ic};
t = num2cell(uniseqCounts);
[t_psmData.psmCountInSample] = t{ic};

%	Create the output struct.
[~, fn, ext] = fileparts(fn_psms);
psmData.meta = struct( ...
	'filename', [fn, ext], ...
	'path', fn_psms, ...
	'params', psmParams ...
	);
psmData.scans = t_psmData;


