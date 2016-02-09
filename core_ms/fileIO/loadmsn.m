function msnData = loadmsn(filesIn)
% LOADMSN load mass spectra from files in the MS1 or MS2 formats.
%	MSNDATA = LOADMSN(FILESIN) loads *.ms1/*.ms2/*.msn file(s) and returns
%	structure array of the imported data.
% 
% MS1, MS2 format spec defined in:
%	 McDonald WH, et al. MS1, MS2, and SQT—three unified, compact, and easily 
%		parsed file formats for the storage of shotgun proteomic spectra and 
%		identifications. Rapid Commun. Mass Spectrom. 2004, 18:2162–2168.
%	http://fields.scripps.edu/sequest/unified/MS2Format.html
% Conversion tools are pretty lax in following the rules, as documented below.
% 
% General structure:
% (Each row is tab-delimited except for space delimited mz-inten pairs.)
%	H - header (required). No spaces/tabs in label.
%	S - scan (required)
%	I - charge-independent vals (optional). No spaces/tabs in label.
%	Z - charge (optional, MS2 only, multiple lines if ambiguous charge vals)
%	D - charge-dependent vals (optional, MS2 only). No spaces/tabs in label.
%	mass/intensity data pairs, space delimited
% 
% MS1 format:
%	H	[label]	[value]
%	S	[first]	[second]
%	I	[label]	[value]
%	[m/z] [intensity]
% 
% MS2 format:
%	H	[label ]	[value ]
%	S	[first ]	[second]	[precursor m/z]
%	I	[label ]	[value ]
%	Z	[charge]	[MH+   ]
%	D	[label ]	[value ]
%	[m/z] [intensity]
% 
% Header fields
%	required: CreationDate, Extractor, ExtractorVersion, ExtractorOptions
%	# These 'required' fields are often not included in ms2 files.
%	optional: IAnalyzer, IAnalyzerVersion, IAnalyzerOptions, InstrumentType,
%			Comment, InstrumentSN
%	optional (MS2): DAnalyzer, DAnalyzerVersion, DAnalyzerOptions, SortedBy
% 
% 
% MSConvert (proteowizard 3.0.4472) doesn't follow formal MS2 file definition.
%	- It puts spaces in the header labels when they should be camel case.
%	- It creates unexpected header field 'Source file'.
%	- It delimits some header fields-val pairs with space, not tab.
%	- It doesn't include required header field 'ExtractorOptions'.
%	- It writes 0 for the precursor m/z value.
% mzxml2msn.jar (1.0.1 release date: 10/26/2009) also fails formal validation.
%	- Doesn't write required header fields: CreationDate, ExtractorOptions
%	- It creates unexpected header fields 'OriginalMzxmlFile', 'MsLevle'.
% Ms2Writer.h (6870 2014-11-03) from MacCoss's BiblioSpec is also out of Spec.
%	- Charge-dependent label 'modified seq' has a space in it.
%	- Doesn't write required header fields: ExtractorVersion, ExtractorOptions
% 
% WARNING: This function assumes the first & second scan numbers are the same.
% @TODO don't make the above assumption.
% @TODO make this faster. (DNmsnScanData.m is faster but more rigid.)
% 
% See also SAVEMSN.

%%

if ~exist('filesIn','var'), filesIn = []; end

% @TODO deprecate io.getfilepath()
[filesIn, fnames] = io.getfilepath(filesIn, '*.ms1;*.ms2;*.ms*', ...
	'select msn file(s)', 'MultiSelect','on');
if isempty(filesIn), return, end

nfiles = numel(filesIn);
msnData(nfiles, 1) = struct;
for n = 1:nfiles
	[msnData(n).meta.fileheader, msnData(n).scans] = load_file( filesIn{n} );
	msnData(n).meta.filename = fnames{n};
	msnData(n).meta.path = filesIn{n};
end

end

%%% Load single MS1/MS2/MSn file.
function [header, scans] = load_file(filepath)
	
	fid = fopen(filepath,'r');
	z_cleanupObj = onCleanup(@() fclose(fid));
	
	%%% read and parse header
	while fscanf(fid,'%c',1) ~= 'S'  % increment fid until first scan line
		fgets(fid); % go to next line
		fidpos = ftell(fid);
	end
	firstScanLine = fgets(fid);
	frewind(fid)
	t = fread(fid, fidpos, '*char')';
	headerstr = regexprep( t, '\r?\n?$','');
	
	%	parse header
	t_split = regexp(regexp(headerstr,'[\r\n]+','split'),'\t','split')';
	nheader = numel(t_split);
	
	header = cell(nheader,2);
	mask = false(nheader,1);
	for i = 1:nheader
		if numel(t_split{i}) > 1 && ...
				strcmp(t_split{i}{1},'H')
			header{i,1} = t_split{i}{2};
			header{i,2} = [t_split{i}{3:end}];
		else
			mask(i) = true;
		end
	end
	header(mask,:) = [];
	
	%	Does the first scan event have a precursor m/z value?
	switch numel(sscanf(firstScanLine,'\t%d'))
		case 2
			msLevel = 1;
		case 3
			msLevel = 2; % 2ndary scans - ms2 or greater
		otherwise
			error('loadmsn:UnknownScanLineFormat', ...
				'Unknown format for scan line in file:  %s',msnData(n).filename)
	end
	
	%%% Read remainder of file
	fseek(fid,2,'cof'); % move past 'S\t' scan line marker.
	fcontent = fread(fid,'*char')'; % @TODO may cause problems if limited memory
	clear z_cleanupObj
	
	% Split raw string into scans & each scan into {S, otherinfo, peaklist} strings
	scanblocks = regexp(fcontent,'\nS','split');
	scanblocks = regexp(scanblocks,'\n','split','once');
	scanblocks = cat(1,scanblocks{:});
	scanblocks(:,2) = regexp(scanblocks(:,2),'\n(?=\d)','split','once');
	scanblocks(:,[2 3]) = cat(1,scanblocks{:,2});
	
	nscans = size(scanblocks,1);
	
	% Interpret the 'S' field strings.
	%	This assumes the first & second scan numbers are always the same.
	if msLevel == 1
		t = num2cell(sscanf([scanblocks{:,1}], '%u\t%*u\t', inf));
		scans = struct('nscan',t);
	else
		t = num2cell(reshape( sscanf([scanblocks{:,1}], '%u\t%*u\t%f\t', inf), 2, []));
		scans = cell2struct(t,{'nscan','pmz'},1);
	end
	
	% Reg expression parses the remainder of each scan header.
	expr_scanHeader = '([IZD])\t(.*?)\t(.*?)(?:$|[\n\r]+)';
	
	% Dynamically build the remaining data for each scan. There's no assumption
	%	that fields are consistent between scans, but they probably should be.
	for i = 1:nscans
		t = regexp(scanblocks{i,2}, expr_scanHeader, 'tokens');
		linedata = cat(1,t{:});
		zinds = 'Z' == [linedata{:,1}];
		if any(zinds)
			t = str2double(linedata(zinds,[2 3]));
			scans(i).pcharge = t(:,1);
			scans(i).pmass = t(:,2);
		end
		idinds = find(~zinds);
		for j = 1:numel(idinds)
			scans(i).(linedata{idinds(j),2}) = linedata{idinds(j),3};
		end
		scans(i).spec = sscanf(scanblocks{i,3},'%f',[2 inf])';
	end
	
end

