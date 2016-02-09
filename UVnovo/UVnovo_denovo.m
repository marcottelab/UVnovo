function varargout = UVnovo_denovo(Meta, varargin)
% UVNOVO_DENOVO
% 
% @TODO add documentation.
% 
% See also UVNOVO, UVNOVO_TRAIN.


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% For development & debugging purposes. ARG1 may be the name of a subfunction
% instead of Meta. In this case, only that subfunction is called and any outputs
% are immediately returned.
if ischar(Meta) && strcmp(which(Meta), which(mfilename))
	fh = str2func(Meta);
	[varargout{1:nargout}] = fh(varargin{:});
	return
end


% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %
% Load spectra and UVnovo trained models.
fprintf(1, 'Importing spectra and UVnovo models ... ')
[msData, Ens, TransMat, filesIn] = import_data(Meta);
fprintf(1, 'Done.\n')


varargout = {};

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Subfunctions


function [msData, Ens, TransMat, filesIn] = import_data(Meta)
	% Load spectra, random forest, and initialize aa mass transition matrix.
	
	paths = Meta.paths;
	filetypes = fieldnames(paths);
	filesIn = cell(3,1);
	if ismember('test', filetypes)
		% Load spectra from test file created by UVnovo_partition, or any
		% serialized MAT file that alreday contains a 'msData' var.
		fn = paths.test.path;
		filesIn{1} = fn;
		s = io.loadSer(fn, 'asStruct');
		if ~isfield(s, 'msData')
			error('UVnovo_denovo:import_data:msDataNotInMatFile', ...
				'File must contain the msData spectra variable.\n\t%s', fn)
		end
		msData = s.msData;
		
	elseif ismember('spectra', filetypes)
		% Import from MS2 file.
		filesIn{1} = paths.spectra.path;
		msData = import_spectra(filesIn{1}, Meta.params.pre.import_spectra);
	end
	
	% Import random forest.
	filesIn{2} = paths.rf.path;
	s = io.loadSer(filesIn{2}, 'asStruct');
	Ens = s.Ens;
	clear s
	
	% Init residue mass transition matrix.
	filesIn{3} = paths.aamodel.path;
	TransMat = massTransMat(filesIn{3}, Meta.params.AAs);
	
end

