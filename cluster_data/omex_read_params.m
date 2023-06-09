function [parout, given] = omex_read_params(params, args)
%[parout, given] = omex_read_params(params, args)
%	Make a parameter struct from variable input parameters.
%
%   Input parameter:
%      args   The varargin cell array of the function call with the names of
%               the parameters and their values in Matlab style, i.e.
%               'Size', 4, 'Method', 'random', ...
%      params The parameters struct containing fields for all possible
%             parameters set beforehand to their default values.
%
%   Output parameter:
%      parout The parameter struct with the named parameters given in args
%             and adjusted to the values given in args plus the default
%               parameters given in params.
%      given  The field names of the actual given parameters in args.
%
%   Comments:
%       If a function uses omex_read_params, but the complete parameter
%       struct is assembled alread (e.g. by the calling programm) you can
%       also give the parameter struct instead of args by calling your
%       functions as:
%           function(usual_parameter, 'params', your_struct)
%
%       If you have the varargin cell array already compiled and just want
%       to pass it to an inner function call your function as:
%           function(usual_parameter, 'varargin', your_varargin);

switch nargin
    case 0
        error('Not enough input arguments!');
    case 1
        args = {};
end

given = {};
parout = params;
size_args = size(args,2);

for i = 1:2:size_args
    if ~ischar(args{i})
        error('at position %d parameter name was a %s, must be a string!', i, class(args{i}));
    end
    if size_args < i+1
        error('no value given for parameter: ''%s''!', args{i});
    end
    if strcmp(args{i}, 'params')
        if ~isstruct(args{i+1})
            error('expected struct as "params" argument')
        end
        flds = fieldnames(args{i + 1});
        for j = 1:size(flds)
            % we do not give any errors (todo)
            parout.(flds{j}) = args{i + 1}.(flds{j});
            given = [given, flds{j}]; %#ok<AGROW>
        end
    elseif strcmp(args{i}, 'varargin')
        parout = omex_read_params(params, args{i + 1});
    elseif isfield(params, args{i})
        parout.(args{i}) = args{i + 1};
        given = [given, args{i}]; %#ok<AGROW>
    else
        error('unknown parameter: ''%s''!', args{i});
    end
end

end
