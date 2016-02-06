
% function takes struct as input and returns struct with field names
%   changed according to fieldnameStandardize.txt
% @TODO Review this code. It's been years since I wrote most of this (<=2012).
function [struct_out, var_names] = fieldnameStandardize(struct_in,no_userprompt)

% If false, prompt user to approve field changes or assign different names.
if ~exist('no_userprompt','var')||isempty(no_userprompt),no_userprompt=1;end
%%
fnames_in = fieldnames(struct_in);
fnames_in(:,2) = {[]};
field_conversions = cell(1000,2);

[pn,fn] = fileparts(mfilename('fullpath'));
fn_conversionTemplate = fullfile(pn, [fn,'.txt']);
if exist(fn_conversionTemplate,'file') ~= 2
		error(exceptionID('FileNotFound',1), 'File %s not found.', fn);
end

fid = fopen(fn_conversionTemplate,'r+');
z_cleanupObj = onCleanup(@() fclose(fid));
c = 1;
while ~feof(fid)    % populate table of field name conversions
    s = fgets(fid);
    if s(1)=='+'
        cur_fname = sscanf(s,'%*c %s',1);
%         field_conversions(c,:) = {cur_fname,cur_fname};
        c = c+1;
    elseif s(1)=='='
        field_conversions(c,:) = {sscanf(s,'%*c %s',1),cur_fname};
        c = c+1;
    end
end
field_conversions(c:end,:) = [];

n_fn = size(fnames_in,1);
for i = 1:n_fn
    field_ind = find(strcmp(fnames_in{i},field_conversions(:,1)));
    if field_ind
        fnames_in(i,2) = field_conversions(field_ind,2);
    end
end
fchange_inds = 1:n_fn;

if ~no_userprompt
    fwidthmax = max(cellfun(@length,fnames_in(:,1)));
    for i = 1:n_fn
        fprintf('%-5.5s %s',[num2str(i) ': '], [fnames_in{i,1} ' ' ...
            repmat('.',1,fwidthmax-length(fnames_in{i,1})+2) ' ' fnames_in{i,2}]);
        fprintf('\n')
    end
    cont=[];
    while isempty(cont)     % choose how to proceed
        t = input('Change all (1), choose set to change (2), or quit (0)? ','s');
        if t=='1'
            cont=1;
        elseif t=='0'
			struct_out = struct_in;
			var_names = fnames_in(:, [1 1]);
            return
        elseif t=='2'
            cont=2;
        end
    end
    
    if cont==2 % choose set of field names to change
        f = regexp(input('choose field names to change (ex. 1 4-8 10 3):  ', 's'),'[,\s]','split');
        fchange_inds = [];
        for k = 1:numel(f)  % parse input into list of files indices
            if ~isnan(str2double(f{k}))
                fchange_inds = [fchange_inds;str2double(f{k})]; %#ok<AGROW>
            elseif numel(regexp(f{k},'-'))==1
                bounds = sort(cellfun(@str2double,regexp(f{k},'-','split')));
                fchange_inds = [fchange_inds;(bounds(1):bounds(2))']; %#ok<AGROW>
            end
        end
        fchange_inds = unique(fchange_inds);
        fchange_inds(~ismember(fchange_inds,1:n_fn)) = [];
    end
end
changemask = false(n_fn,1);
changemask(fchange_inds) = true;
ind_newnames = changemask & ~cellfun('isempty',fnames_in(:,2));

if any(ind_newnames)
	fnames_out = fnames_in(:,1);
	fnames_out(ind_newnames) = fnames_in(ind_newnames,2);
	
	c_struct_in = struct2cell(struct_in);
	struct_out = cell2struct(c_struct_in, fnames_out, 1);
	var_names = [fnames_in(:,1), fnames_out];
else
	struct_out = struct_in;
	var_names = fnames_in(:,[1 1]);
end

