function paramsOut = initParams(paramsDefault, paramsUser, varargin)
% INITPARAMS Merge default and user-defined parameters.
%	Combine a struct of default parameters with user-defined
%		parameters, with the later taking precedence. The user-defined params 
%		can be a struct, or name/value pairs which may be encapsulated in a
%		cell.
% 
%	INITPARAMS(PARAMSDEFAULT, PARAMSUSER, VARARGIN)
%	INITPARAMS(PARAMSDEFAULT, VARARGIN)
% 
% INPUT ARGS
%	paramsDefault: <struct> default parameters
%	paramsUser: <struct OR cell> user provided parameters.
%		If <struct>, must have same numel as default when both are non-scalar.
%		A paramsUser <cell> can be 1-D or {n x 2} of {param,val; ...}.
%	varargin: name/value pairs of additional user parameters.
% 
% OUTPUT
%	paramsOut: <struct> structured like paramsDefault with additional fields or
%		dimensions depending on paramsUser and/or varargin.
% 
% EXAMPLES
%		paramsDef = struct('p1',1, 'p2',{{'hallo','wrold'}});
%		s1 = struct('p1',2, 'p0','yes', 'p55',[2,3,55]);
%		s2 = struct('p1',2, 'p0','NO');
%		c1 = {'p1', 2, 'p0', 'yes', 'p55', [2,3,55]};
%	Each of these is equivalent:
%		p0 = initParams(paramsDef, s1);
%		p1 = initParams(paramsDef, c1);
%		p2 = initParams(paramsDef, 'p1', 2, 'p0', 'yes', 'p55', [2,3,55]);
%		p3 = initParams(paramsDef, s2, 'p0', 'yes', 'p55', [2,3,55]);
%		p4 = initParams(paramsDef, {'p1', 2, 'p0', 'yes'}, 'p55', [2,3,55]);
%		isequal(p0, p1, p2, p3, p4)
% 

if ~isstruct(paramsDefault)
	error('initParams:defParamsNotStruct',...
		'First argument must be struct specifying default parameters.')
end

if ~exist('paramsUser','var') || ...
		( isempty(paramsUser) && numel(varargin)==0 )
	paramsOut = paramsDefault;
	return
end

if ~(isstruct(paramsUser) || iscell(paramsUser) || ...
		(ischar(paramsUser) && numel(varargin)) )
	error('initParams:invalidUserParams', ...
		'Argument 2 must be empty or specify user params.')
end


if iscell(paramsUser)
	% Put user params into a struct.
	paramsUser = paramsUser';
	paramsUser = cell2struct(paramsUser(2:2:end), paramsUser(1:2:end), 2);
end

if isstruct(paramsUser) && numel(varargin)
	% Merge ARG_2 user params & name/val pairs. Name/val pairs take precedence.
	paramsUser = initParams(paramsUser,varargin);
end

if ischar(paramsUser)
	% Convert ARG_2 and varargin name/value pairs into struct.
	pairs = [paramsUser; varargin(:)]';
	try
		paramsUser = cell2struct(pairs(2:2:end), pairs(1:2:end), 2);
	catch me
		error('io:initParams:invalidNameValArgs', me.message)
	end
end


% Merge default and user parameters.
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
		assert(nDef == nUser, 'io:initParams:DefAndUserNumelMismatch', ...
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
				t_params(o).(field) = initParams( ...
					t_params(o).(field), ...
					paramsUser(iu).(field) );
			else
				t_params(o).(field) = paramsUser(iu).(field);
			end
		end
	end
end

paramsOut = t_params;
