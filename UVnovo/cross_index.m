function mapi = cross_index(psmData, msData, mapi)
% CROSS_INDEX index between data structures.
%	This function is very specific and only maps between psmData.scans and
%	msData.scans. It would be more useful if generalized to index between
%	arbitrary data structures.
% 
% INPUT ARGUMENTS
%	psmData
%	msData
%	mapi: <struct> Optional. If provided, this function updatesa with additional
%		mappings. !! This is 
% OUTPUT
%	mapi: <struct> index mappings between the provided data structures.
% 
% 
% 
% @TODO Generalize and abstract out any specific data structure needs.
% @TODO testing & validation. Make sure all psms & scans are matched.


% ismember() is fast enough to just call twice, rather than writing more
%	efficient code. That's easy to change if this function ever gets fleshed out
%	& generalized.
[~, psm2ms] = ismember([psmData.scans.nscan]', [msData.scans.nscan]');
[~, ms2psm] = ismember([msData.scans.nscan]', [psmData.scans.nscan]');

mapi.psmData.msData = psm2ms;
mapi.msData.psmData = ms2psm;

