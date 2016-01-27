function paramsOut = init_params(paramsDefault,paramsUser,varargin)
% INIT_PARAMS Merge default and user-defined parameters.
%	Combine programatically-defined default parameters with user-defined
%		parameters, with the later taking precedence. The user-defined params 
%		can be a struct, name/value pairs, or a combination of both. The 
%		name/value pairs may be encapsulated in a cell or not.
%	
% INPUT ARGS
%	paramsDefault: <struct> default parameters
%	paramsUser: <struct, cell, OR str> user provided parameters
%		If <struct> must have same numel as default when both are non-scalar.
%		<cell> can be 1-D or {n x 2} of {param,val; ...}
%	varargin: name/value pairs (paramsUser will contain name for first pair).
% 
% OUTPUT
%	paramsOut: <struct> structured like paramsDefault with additional fields or
%		dimensions depending on paramsUser and/or varargin.

if ~exist('paramsUser','var'), paramsUser = []; end

if ~isstruct(paramsDefault)
	error('init_params:defParamsNotStruct',...
		'First argument must be struct specifying default parameters.')
end

if isempty(paramsUser) && numel(varargin) == 0
	paramsOut = paramsDefault;
	return
end

if ~(isstruct(paramsUser) || iscell(paramsUser) || ...
		(ischar(paramsUser) && numel(varargin)) )
	error('init_params:invalidUserParams', ...
		'Argument 2 must be empty or specify user params.')
end


%%%	Put all user params into a struct if they aren't already.
if isstruct(paramsUser) && numel(varargin)
	% Merge user struct & name/val pairs. Name/val pairs take precedence.
	paramsUser = init_params(paramsUser,varargin);
elseif numel(varargin) || iscell(paramsUser)
	% Convert name/value argument pairs into struct.
	paramsUser = paramsUser';
	pairs = [paramsUser(:); varargin]';
	try
		paramsUser = cell2struct(pairs(2:2:end), pairs(1:2:end), 2);
	catch me
		error( struct(  'identifier','io:init_params:invalidNameValArgs', ...
						'message',me.message ))
	end
end

%%% Merge default and user parameters.
nDef = numel(paramsDefault);
nUser = numel(paramsUser);
nOut = max(nUser, nDef);

switch (nDef>1) + 2*(nUser>1)
	case {0,1} % both == 1 OR nDef > 1
		t_params = paramsDefault;
		indsUser = ones(nOut,1);
	case 2 % nUser > 1
		t_params = repmat(paramsDefault,size(paramsUser));
		indsUser = 1:nOut;
	case 3 % both > 1
		assert(nDef == nUser, 'io:init_params:DefAndUserNumelMismatch', ...
			'Default and user parameter structs must have == numel when both not scalar.')
		t_params = paramsDefault;
		indsUser = 1:nOut;
end

p_fields = fieldnames(paramsUser);
for i = 1:numel(p_fields)
	field = p_fields{i};
	if ~isfield(t_params, field)
		[t_params.(field)] = paramsUser(indsUser).(field);
	else % update existing field with user-defined value
		for o = 1:nOut
			iu = indsUser(o);
			if isstruct(paramsUser(iu).(field)) && isstruct(t_params(o).(field))
				% recursively fill in parameters of nested structs
				t_params(o).(field) = init_params( ...
					t_params(o).(field), ...
					paramsUser(iu).(field) );
			else
				t_params(o).(field) = paramsUser(iu).(field);
			end
		end
	end
end

paramsOut = t_params;
