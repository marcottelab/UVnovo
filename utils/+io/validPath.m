function [pathOut, attribs] = validPath(pathIn)
% VALIDPATH get longest accessible path from pathIn.
%	If pathIn is not valid, shortens by a single level until it's valid.
%	This returns pwd if root of pathIn doesn't exist.

t_path = pathIn;
t_len0 = length(t_path);
while ~exist(t_path,'dir')
	t_path = fullfile(t_path,['..' filesep]);
	t_len = length(t_path);
	if t_len > t_len0
		t_path = pwd; break
	end
	t_len0 = t_len;
end

[~,attribs] = fileattrib(t_path);
pathOut = attribs.Name;
