function scans = MSfilterScans(scans, paramsIn)
% MSfilterScans: Apply basic filtering to mass spectra.
%	scans = MSfilterScans(scans, paramsIn)
% ARGS
%	scans: <cell array of [m x 2]> {[peak m/z, peak intensity];...}
%	paramsIn: (optional) <struct>
% OUTPUT
%	scans: <cell array of [m x 2]> filtered spectra
% 
% This takes a set of scans & removes peaks not meeting user-defined criteria.
% Named filters can be added as needed.

paramsDef = struct( ...
	'minInten', [], ... % remove peaks below this intensity
	'removeZeroInten', true, ... % remove peaks with zero intensity
	'massRange', [] ... % mass range to return. NOT IMPLEMENTED!!!
	);
if ~exist('paramsIn','var'), paramsIn = []; end
params = init_params(paramsDef, paramsIn);

for i = 1:numel(scans);
	
	if ~isempty(params.minInten) && params.minInten > 0
		scans{i}( scans{i}(:,2) < params.minInten, : ) = [];
		
	elseif params.removeZeroInten
		scans{i}( scans{i}(:,2) == 0, : ) = [];
		
	end
	
end

