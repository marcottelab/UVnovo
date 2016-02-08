function paramsOut = initParams(paramsDefault, paramsUser, varargin)
% INITPARAMS update a defaults parameter struct with user-defined values.
%	Combine a struct of default parameters with user-defined parameters, with
%	the later taking precedence. The user-defined params can be a struct or
%	name/value pairs. The pairs may be encapsulated in a cell array. Additional
%	params may be passed in varargin as name/val pairs, these taking precedence.
% 
%	INITPARAMS(PARAMSDEFAULT, PARAMSUSER, VARARGIN)
%	INITPARAMS(PARAMSDEFAULT, VARARGIN)
% 
% INPUT ARGS
%	paramsDefault: <struct> default parameters
%	paramsUser: <struct OR cell> user provided parameters.
%		A paramsUser <cell> can be 1-D or {n x 2} of {param,val; ...}.
%	varargin: name/value pairs of additional user parameters.
% 
% OUTPUT
%	paramsOut: <struct [size(paramsUser)]> fields structured like paramsDefault
%		with additional fields depending on paramsUser and/or varargin.
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
% 
%	Default and user params needn't be scalars. The returned parameters struct
%	takes the same size as paramsUser. When paramsDefault is larger, elements
%	with a greater index than nUser are discarded.
% 
%	NOTE: When paramsDefault is a nonscalar struct with fewer elements than
%		paramsUser, paramsDefault is padded to be the same size, using the last
%		element of paramsDefault to pad with.
%	Also, be careful if user and default param structs are both multidimensional
%		with different sizes. The individual elements may not align correctly.


if ~isstruct(paramsDefault)
	error('initParams:defParamsNotStruct',...
		'First argument must be struct specifying default parameters.')
end

if ~exist('paramsUser','var') || ...
		( isempty(paramsUser) && numel(varargin)==0 )
	paramsOut = paramsDefault;
	return
end

if isempty(paramsUser)
	paramsUser = varargin;
	varargin = {};
end

if ~(isstruct(paramsUser) || iscell(paramsUser) || ...
		(ischar(paramsUser) && numel(varargin)) )
	error('initParams:invalidUserParams', ...
		'Argument 2 must be empty or specify user params.')
end


if iscell(paramsUser)
	% Put user params into a struct.
	if size(paramsUser,2) == 2
		paramsUser = paramsUser';
	end
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
nOut = nUser;

% Initilize temp struct of defaults that will be updated with user params.
if nDef >= nUser
	t_params = reshape(paramsDefault(1:nUser), size(paramsUser));
else
	t_params = repmat(paramsDefault(end), size(paramsUser));
	t_params(1:nDef) = paramsDefault;
end
indsUser = 1:nOut;


p_fields = fieldnames(paramsUser);
for i = 1:numel(p_fields)
	field = p_fields{i};
	if ~isfield(t_params, field)
		[t_params.(field)] = paramsUser(indsUser).(field);
	else % Update existing field with user-defined value.
		for o = 1:nOut
			iu = indsUser(o);
			if isstruct(paramsUser(iu).(field)) && isstruct(t_params(o).(field))
				% Recursively fill in parameters of nested structs.
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
