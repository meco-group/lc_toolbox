function [ lambda ] = buildKnotSeq(varargin)

required = {'breakpoints','degree'};
optional = {};

calls = varargin(1:2:end);
assert(mod(nargin,2) == 0,'Odd number of arguments.');
assert(iscellstr(calls),'The names of your name/value pairs should be strings.');
assert(~any(~ismember(required,calls)),['These arguments are required: ' strjoin(required)]);
for i = 1:nargin/2
    if strcmp(varargin{2*i-1},'breakpoints'); assert(~exist('breakpoints','var'),'You cannot define breakpoints twice.'); breakpoints = varargin{2*i}; end;
    if strcmp(varargin{2*i-1},'degree'); assert(~exist('degree','var'),'You cannot define degree twice.'); degree = varargin{2*i}; end;
end

lambda = [min(breakpoints)*ones(1, degree), sort(breakpoints), max(breakpoints)*ones(1, degree)]; 

end
