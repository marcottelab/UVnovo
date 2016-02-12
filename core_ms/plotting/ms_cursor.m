function output_txt = ms_cursor(~,event_obj)
% Display the position of the data cursor
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).
% 
% Example usage:
% 	h = figure;
% 	dcm_obj = datacursormode(h);
% 	set(dcm_obj,'UpdateFcn',@ms_cursor)


pos = get(event_obj,'Position');
output_txt = {['m/z: ',num2str(pos(1),7)] ['inten: ',num2str(pos(2),7)]};

% If there is a Z-coordinate in the position, display it as well.
if length(pos) > 2
    output_txt{end+1} = ['Z: ',num2str(pos(3),4)];
end

