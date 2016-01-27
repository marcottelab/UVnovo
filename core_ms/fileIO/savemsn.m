function savemsn(fnOut, spec, scanInfo, expHeader, paramsIn)
% SAVEMSN save data in MS1 or MS2 format
%	See loadmsn.m for specification on proper MS1/MS2 format.
% 
% fnOut: file name. Prompt appears if not provided.
% spec: <cell array [n x 1]> n spectra, each a [m x 2] mass/inten array
% scanInfo: <struct array> elements correspond to spec cells
%	REQUIRED FIELDS:
%		scan number (default field name: 'nscan')
%	REQUIRED FOR MS2:
%		precursor m/z (default: 'pmz')
%		precursor charge (default: 'pcharge')
%		precursor mass (default: 'pmass')
%	Provide translator in paramsIn if non-default field names (template below).
% expHeader: <scalar struct>
% paramsIn: <struct>
%	specialFields: <struct> define names for required scanInfo fields.
%	printFields: <cellstr|logical> list of scanInfo fields to write to msn file.
%		<empty(default) | true> print all non-required scanInfo fields.
%		<false> don't print any non-required fields.
%	scanPrecision: <int> precision of sprintf conversion specifier.
%		Since scan data is typically single-precision, any digits beyond 7 or 8 
%		are artifacts of converting from single -> double precision.
% 
% Compared to full ms1/ms2 file specifications:
%	this does not verify that 'required' header fields are provided.
%	this does not yet allow for multiple charge states for a single scan.
%	this treats all scan data as charge-independent (I). No (D) fields written.
% 
% See also LOADMSN.

%%
if ~exist('paramsIn','var') || isempty(paramsIn), paramsIn = []; end
paramsDef = struct(...
	'specialFields', struct(... % translators for required field names
	...  Required field     Fieldname in 'scanInfo'
		'scan_num',        'nscan',...
		'precursor_mz',    'pmz',...
		'precursor_charge','pcharge',...
		'precursor_mass',  'pmass'...
		), ...
	'printFields', {''},... % specify fields to print for each scan
	'scanPrecision', 8 ...  % beyond ~8 is artifact of single -> dbl conversion
	);
params = init_params(paramsDef, paramsIn);
specials = params.specialFields;
frmtprecision = sprintf('%d',params.scanPrecision);

nspec = numel(spec);
assert(nspec == numel(scanInfo), 'savemsn:SpecDataMismatch', ...
	'number of spectra must match number of scanInfo elements')

if ~exist('fnOut','var') || isempty(fnOut)
	if isfield(scanInfo, {specials.precursor_mz})
		fnSpec = {'*.ms2';'*.ms1';'*.*'};
	else
		fnSpec = {'*.ms1';'*.ms2';'*.*'};
	end
	[fn,pn] = uiputfile(fnSpec,'Save file');
	if ~fn, return, end
	fnOut = fullfile(pn,fn);
end

%%
t_dataH = struct2cell(expHeader);
t_dataH(cellfun('isempty',t_dataH)) = {''};
t_indConvert = cellfun(@isnumeric,t_dataH);
t_dataH(t_indConvert) = sprintfc('%.15g',[t_dataH{t_indConvert}]);
dataH = [fieldnames(expHeader), t_dataH]';

frmtH = 'H\t%s\t%s\n';

scanNums = {scanInfo.(specials.scan_num)};
if isfield(scanInfo,specials.precursor_mz)
	scanPmz = {scanInfo.(specials.precursor_mz)};
	frmtS = 'S\t%d\t%d\t%.15g\n';
	dataS = [scanNums;scanNums;scanPmz];
elseif any(isfield(scanInfo, {specials.precursor_charge, specials.precursor_mass}))
	frmtS = 'S\t%d\t%d\t0\n';
	dataS = [scanNums;scanNums];
else
	frmtS = 'S\t%d\t%d\n';
	dataS = [scanNums;scanNums];
end

if islogical(params.printFields) && ~params.printFields
	frmtI = '';
	dataI = {};
else
	fieldsSD = fieldnames(scanInfo);
	
	if ~isempty(params.printFields) && ~islogical(params.printFields)
		printFields = params.printFields;
	else % print all but special fields
		printFields = fieldsSD(~ismember(fieldsSD, struct2cell(specials)));
	end
	
	[~,indI] = ismember(printFields,fieldsSD);
	fieldsI = fieldsSD(indI);
	nI = numel(indI);
	dataI = cell(nI,nspec);
	for f = 1:nI
		dataI(f,:) = {scanInfo.(fieldsI{f})};
	end
	
	t_indConvert = cellfun(@isnumeric,dataI);
	dataI(t_indConvert) = sprintfc('%.15g',[dataI{t_indConvert}]);
	t_frmtI = cell(nI,1);
	for f = 1:nI
		t_frmtI{f} = ['I\t' fieldsI{f} '\t%s\n'];
	end
	frmtI = [t_frmtI{:}];
end

if isfield(scanInfo, specials.precursor_charge)
	pcharge = {scanInfo.(specials.precursor_charge)};
	if isfield(scanInfo, specials.precursor_mass)
		pmass = {scanInfo.(specials.precursor_mass)};
	elseif isfield(scanInfo, specials.precursor_mz)
		pcharge = {scanInfo.(specials.precursor_charge)};
		try
			pmz = [scanInfo.(specials.precursor_mz)];
			assert(numel(pmz) == nspec)
			t = [pcharge{:}];
			pmass = num2cell(pmz.*t-(t-1)*CONSTS.mProton);
		catch
			pmass = cell(1,nspec);
			pmz = {scanInfo.(specials.precursor_mz)};
			indNE = ~cellfun('isempty',pmz) & ~cellfun('isempty', pcharge);
			indNums = cellfun(@isnumeric,pmz(indNE)) & ...
					  cellfun(@isnumeric,pcharge(indNE));
			indNE(indNE) = indNums;
			pmass(indNE) = num2cell(pmz2mh([pmz{indNE}], [pcharge{indNE}])');
		end
	else
		pmass = zeros(size(pcharge));
	end
	dataZ = cat(1,pcharge,pmass);
	frmtZ = ['Z\t%d\t%.' frmtprecision 'g\n'];
else
	frmtZ = '';
	dataZ = {};
end

%	S	[first]	[second]	[precursor m/z]
%	I	[label]	[value]
%	Z	[charge]	[MH+]
%	[m/z] [intensity]

frmtSIZ = [frmtS, frmtI, frmtZ];
dataSIZ = cat(1,dataS,dataI,dataZ);

frmtScan = ['%.' frmtprecision 'g %.' frmtprecision 'g\n'];
%%
fid = fopen(fnOut,'W');
fprintf(fid,frmtH,dataH{:});
for i = 1:nspec
	fprintf(fid,frmtSIZ,dataSIZ{:,i});
	fprintf(fid,frmtScan,spec{i}');
end
fclose(fid);

