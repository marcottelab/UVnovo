function output_txt = UserData_cursor(~,event_obj)
% USERDATA_CURSOR cursor mode for custom UserData tooltips.
%	This one allows non-numeric tip data.
% 
% EVENT_OBJ    Handle to event object - internal var Matlab takes care of
% OUTPUT_TXT   Data cursor text string (string or cell array of strings).
% 
% Usage example:
%	data_object_handle = plot( ... );
%	set(data_object_handle, ...
%		'UserData', struct( 'data',   {Data}, ...
%							'format', {DataFormat}
%			)
%		);
%	dcm_obj = datacursormode(fig_handle);
%	set(dcm_obj,'UpdateFcn',@UserData_cursor)
% Vars in above example:
%	Data: <cell {m}> m: number of vars to show in datatip
%						Each cell contains array with same size as graph data.
%	DataFormat: <cell (1 x m)> m: sprintf-style format for each datatip var
% 
% see also DATACURSORMODE.
% 
% Andrew Horton, 2014

% Get the supplied data struct
userData = get(event_obj.Target,'UserData');
data = userData.data;
dataFormat = userData.format;

% Get index of selected point
try
	sizeData = size(get(event_obj.Target,'CData'));
catch
	sizeData = size(get(event_obj.Target,'XData'));
end
dataSubs = fliplr(get(event_obj,'DataIndex'));
if numel(dataSubs) == 2
	dataInd = sub2ind( sizeData, dataSubs(1), dataSubs(2) );
else
	dataInd = dataSubs;
end

nTags = numel(dataFormat);
output_txt = cell(1,nTags);
for i = 1:nTags
	if iscell(data{i})
		output_txt{i} = sprintf(dataFormat{i}, data{i}{dataInd});
	else
		output_txt{i} = sprintf(dataFormat{i}, data{i}(dataInd));
	end
end	
