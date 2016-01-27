function [fullfiles, fnames] = getfilepath(filesIn, varargin)
% getfilepath: check that input file paths exist or prompt if none provided.
%	varargin: args that may be passed to uigetfile()

if ~exist('filesIn','var') || isempty(filesIn)
	[fnames, pname, filtindex] = uigetfile(varargin{:});
	if filtindex == 0
		fullfiles = {};
		warning('io:getfilepath:NoFilesChosen','No file chosen.')
		return
	end
	if ~iscell(fnames)
		fnames = {fnames};
	else
		fnames = fnames(:);
	end
	nfiles = numel(fnames);
	fullfiles = cell(nfiles,1);
	for i = 1:nfiles
		fullfiles{i} = fullfile(pname, fnames{i});
	end
else
	if ~iscell(filesIn)
		filesIn = {filesIn};
	end
	nfiles = numel(filesIn);
	fullfiles = cell(nfiles,1);
	for i = 1:nfiles
		fn = filesIn{i};
		if ~ischar(fn)
			error(exceptionID('FileNameMustBeString',1), ...
				'File names must be strings.');
		end
		if exist(fn,'file') ~= 2
			error(exceptionID('FileNotFound',1), ...
				'File %s not found.', fn);
		else
			[~,t] = fileattrib(fn);
			fullfiles{i} = t.Name;
		end
	end
	
	if nargout > 1
		fnames = cell(nfiles,1);
		for i = 1:nfiles
			[~, fn, fnExt] = fileparts(fullfiles{i});
			fnames{i} = [fn,fnExt];
		end
	end
end


