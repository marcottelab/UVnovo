function [spepdata, pFields, fnames] = load_xlsx_psms(fn_psm_xlsx, paramsIn)
% LOAD_XLSX_PSMS load Thermo Proteome Discoverer PSMs from Excel spreadsheets.
% 
% INPUT
%	fn_psm_xlsx: filename(s) of psm files to load (*.XLS or *.XLSX format).
%		Files should be formatted like Proteome Discoverer, where each row
%		contains various information relating to a single PSM. Multiple PSM rows
%		can correspond with a single spectrum.
%	paramsIn: <struct> see section below for more parameter documentation.
%		.pFields: <n x 3 cell> define how to interpret columns in psm files.
% 
% OUTPUT
%	spepdata: 
%	pFields: list of fields in spepdata
%	fnames
%
% This function is a vestige of earlier days and is a complete mess!!!
% @TODO Tidy this up significantly or replace.
% 
% See also IMPORT_PSMS.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------OPTIONS------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
paramsDef = struct(...
	'trim_accessions',	0,...	% only keep the first accession listed
	'rank_limit',		3,...	% maximum rank to include
	'add_ppm_adj',		0,...   % add field 'ppm_adj' which is ppm minus avg ppm of Unambiguous PSMs without 'x' in sequence
	'decoy_flag',		[],...	% ex. {'shuffled_';'full-shuff_conJ_1_'};
	'verbosity',        1 ...	% [0-2] Print info messages? 0:none to 2:all.
	);

%%% Specify columns to import, what to name them in the return struct, and which 
%%% of potentially many values for a PSMs to retain. 
%	HEADER: column header name after calling Matlab's genvarname(). The name is
%		not case sensitive. For fields that may take one of several names, each
%		is included in a cellstr. The first match in the cellstr is used.
%	NAME: field name in output psmData struct
%	SELECTION: psms are collapsed based on scan. When there are multiple psms,
%		either (0) all values or (1) the value of the top-ranked PSM, may be
%		returned for each scan. Option 0 returns a cell arrays.
paramsDef.pFields = {
	%HEADER                      NAME           SELECTION
	'FirstScan'                 'scan'          1;...
	'Sequence'                  'sequence'      0;...
	'XCorr'                     'XCorr'         0;...
	{'DeltaScore','deltaScore','x0x0394Score'}		'dScore'	1;...
	{'DeltaCn','deltaCn','x0x0394Cn'}				'dCn'		0;...
	{'DeltaMass0x5BPPM0x5D','x0x0394M0x5Bppm0x5D'}	'ppm'		0;...
	'Rank'                      'rank'          0;...
    'Modifications'             'PTMs'          0;...
	'm0x2Fz0x5BDa0x5D'          'pmz'           1;...
	'MH0x2B0x5BDa0x5D'          'pmass'         1;...
	'Charge'                    'pcharge'       1;...
	'RT0x5Bmin0x5D'             'RT'            1;...
	'ConfidenceLevel'           'confidence'    0;...
% 	'PSMAmbiguity'              'ambiguity'     0;...
% 	'x0x23Proteins'             'nProteins'     0;...
% 	'x0x23ProteinGroups'        'nProtGroup'    0;...
% 	'ProteinGroupAccessions'    'accession'     0;...
% 	'SearchEngineRank'          'SearchEngineRank'  0;...
	'q0x2DValue'                'qValue'        0;...
	'PEP'                       'PEP'           0;...
% 	'IsolationInterference0x5B0x250x5D'     'isolationInterference'   1;...
% 	'IonInjectTime0x5Bms0x5D'   'InjectTime'    1;...
% 	'Intensity'                 'Intensity'     1;...
%	'PeptidesMatched'           'pepsMatched'   0;...
% 	'Probability'               'probability'   0;...
% 	'x0x23MissedCleavages'      'nMissedCleavages'  0;...
% 	'LastScan'                  'lastScan'      1;...
%	'SpectrumFile'              'spectrumFile'  1;...
% 	'ActivationType'            'activation';...
% 	'MSOrder'                   'MS_Order';...
% 	'IonsMatched'               'ionsMatched'   0;...
% 	'MatchedIons'               'matchedIons'   0;...
% 	'TotalIons'					'totalIons'     0;...
	};

if ~exist('paramsIn','var'), paramsIn = []; end
params = init_params(paramsDef, paramsIn);

verbosity = params.verbosity;
pFields = params.pFields;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------MAIN-------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('fn_psm_xlsx','var') || isempty(fn_psm_xlsx)
	[fnames, pn] = uigetfile('*.xlsx','MultiSelect','on');
	try if ~fnames, disp('no files chosen'), return, end, catch, end
	if ~iscell(fnames)
		fnames = {fnames};  %convert file name to cell if it isn't
	end
	fpaths = cellfun(@(fnames)fullfile(pn,fnames),fnames,'UniformOutput',false);
else
	if iscell(fn_psm_xlsx)
		fpaths = fn_psm_xlsx;
	else
		fpaths = {fn_psm_xlsx};
	end
end

n_files = length(fpaths);     %number of files
spepdata(n_files).sample = 0;
for n = 1:n_files
	[~, fn_bare, ext] = fileparts(fpaths{n});
	fn = [fn_bare,ext];
	spepdata(n).fname = fn;
	
	printf_msg( 1, 'Loading file: %s\n', fn)
	
	tic_xlsread = tic;
	[~, ~, pepdata] = xlsread(fpaths{n},1);
	
	printf_msg( 2, '\tFile loaded in %g seconds.\n\tParsing...\n', ...
		toc(tic_xlsread))
	
	tic_parse = tic;
	clear uselessvariable
	headers = pepdata(1,:);
	pepdata(1,:) = [];
	nan_headers = ~cellfun(@(x) any(~isnan(x)),headers);
	headers(nan_headers) = [];
	pepdata(:,nan_headers) = [];
	
	field_inds = zeros(size(pFields,1),1);
	nfields = size(pFields,1);
	headers_normalized = lower( genvarname(headers) );
	for i = 1:size(pFields,1)
		if ~iscellstr(pFields{i}) || numel(pFields{i}) == 1
			% Fields with only 1 possibility in headers.
			t_hit = strcmpi(pFields{i},headers_normalized);
			
		else
			% Fields with multiple potential column names.
			%	In case of multiple pField<==>header matches, the first matching
			%	pFields{i} element is selected. pFields{i} is changed to include
			%	just that name.
			j = 1;
			t_hit = [];
			while j<=numel(pFields{i}) && ~any(t_hit)
				t_hit = strcmpi(pFields{i}{j},headers_normalized);
				j = j+1;
			end
			if any(t_hit)
				pFields{i} = pFields{i}{j-1};
			end
		end
		
		if nnz(t_hit)
			field_inds(i) = find(t_hit, 1);
		end
	end
	
	% @TODO 'finding' the field indices this way is gross. Clean it up!
	ind_rank = find(strcmpi('Rank',headers_normalized));
% 	ind_scan = find(strcmpi('FirstScan',varname_headers));
	[~, ia] = ismember(lower(pFields(:,2)),{'scan','nscan'});% HACK, all of this
	ind_scan = find(strcmpi(pFields{ia == min(ia(ia>0)),1},headers_normalized));
	ind_acc = find(strcmpi('ProteinGroupAccessions',headers_normalized));
	inds_0 = [pFields{:,3}]==0 & logical(field_inds)';     % indices of fields containing all values
	inds_1 = [pFields{:,3}]==1 & logical(field_inds)';     % indices of fields containing just first value
	pepdata(isnan([pepdata{:,ind_scan}]),:) = [];   % remove rows that don't have a scan number (most likely empty)
	
	if params.trim_accessions
		try
			pepdata(:,ind_acc) = regexprep(pepdata(:,ind_acc),';[\S\s]*','');
		catch
			pepdata(~cellfun(@ischar,pepdata(:,ind_acc)),ind_acc) = {'_x_'};
			pepdata(:,ind_acc) = regexprep(pepdata(:,ind_acc),';[\S\s]*','');
		end
	end
	if exist('params.decoy_flag','var') && ~isempty(params.decoy_flag)
		pepdata = [pepdata num2cell(zeros(size(pepdata,1),1))];     %#ok<AGROW>
		pFields = [pFields; {'not_in_excel' 'decoy' 0}];            %#ok<AGROW>
		inds_0 = [inds_0 true];                                     %#ok<AGROW>
		field_inds = [field_inds; size(pepdata,2)];                 %#ok<AGROW>
		for i = 1:numel(params.decoy_flag)
			inds_decoy = strncmpi(params.decoy_flag{i},pepdata(:,ind_acc),size(params.decoy_flag{i},2));
			pepdata(:,end) = num2cell([pepdata{:,end}]' | inds_decoy);
			pepdata(:,ind_acc) = regexprep(pepdata(:,ind_acc),params.decoy_flag{i},'d_');
		end
	end
	
	pepdata = sortrows(pepdata,ind_rank);
	pepdata([pepdata{:,ind_rank}]'>params.rank_limit,:) = []; % deletes all with rank>rank_limit
	pepdata = sortrows(pepdata,ind_scan);
	
	scannums = unique([pepdata{:,ind_scan}]);
	nscans = numel(scannums);
	
	j = 1;
	% prevent j from exceeding pepdata dimension
	pepdata = [pepdata; cell(1,size(pepdata,2))];   %#ok<AGROW>
	o = ones(1,nnz(inds_0));
	scandata = cell(nfields,nscans);
	% spepdata(nscans).scan = 0;
	for i = 1:nscans    % compress pepdata for each scan
		j_start = j;
		while pepdata{j,ind_scan}==scannums(i), j = j+1; end
		%     s(i).seq = {pepdata{j_start:j-1,2}}';
		scanfull = pepdata(j_start:j-1,field_inds(inds_0));
		scanfirst = pepdata(j_start,field_inds(inds_1));
		%     scanrest = pepdata(j_start+1:j-1,field_inds(inds_0)); % for future use if needed
		scandata(inds_0,i) = mat2cell(scanfull,size(scanfull,1),o);
		scandata(inds_1,i) = scanfirst;
	end
	structpepdata = cell2struct(scandata,pFields(:,2));
	clear scandata
	
	if params.add_ppm_adj	% add field ppm_adj (ppm adjusted by subtracting median)
		seqs = cat(1,structpepdata.sequence);
		ind_notx = cellfun(@isempty,regexp(seqs(:,1),'x'));
		ppm = cat(1,structpepdata.ppm);
		if isfield(structpepdata,'ambiguity')
			ambiglist = cat(1,structpepdata.ambiguity);
			indunambig = strcmp('Unambiguous',ambiglist);
			medianppm = median([ppm{indunambig&ind_notx}]);
		else
			medianppm = median([ppm{ind_notx}]);
		end
		for j = 1:numel(structpepdata)
			structpepdata(j).ppm_adj = num2cell(cell2mat(structpepdata(j).ppm)-medianppm);
		end
	end
	spepdata(n).sample = structpepdata;
	clear structpepdata seqs ppm
	
	printf_msg(2, '\tParsed in %g seconds\n', toc(tic_parse))
end

if params.add_ppm_adj
	pFields(end+1,:) = {'...internalCalc...','ppm_adj',0};
end

fnames = fpaths; % to be consistent with old version. @.DEP

printf_msg(2,'ALL DONE\n')

%%%%%%%%%%%%%%%%%%%%
%%% Helper functions

	function printf_msg(msg_importance, format_str, varargin)
		% Format and print string to command window if sufficient verbosity.
		if verbosity >= msg_importance
			fprintf(1, format_str, varargin{:})
		end
	end

end

