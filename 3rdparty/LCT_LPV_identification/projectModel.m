function [ regMdlProj ] = projectAnalyticModel( varargin )

required = {'model','schParam','slackVars','threshold'};
optional = {};

calls = varargin(1:2:end);
assert(mod(nargin,2) == 0,'Odd number of arguments.');
assert(iscellstr(calls),'The names of your name/value pairs should be strings.');
assert(~any(~ismember(required,calls)),['These arguments are required: ' strjoin(required)]);
for i = 1:nargin/2
    if strcmp(varargin{2*i-1},'model'); assert(~exist('model','var'),'You cannot define model twice.'); model = varargin{2*i}; end;
    if strcmp(varargin{2*i-1},'schParam'); assert(~exist('schParam','var'),'You cannot define schParam twice.'); schParam = varargin{2*i}; end;
    if strcmp(varargin{2*i-1},'slackVars'); assert(~exist('slackVars','var'),'You cannot define griddedFRFs twice.'); slackVars = varargin{2*i}; end;
    if strcmp(varargin{2*i-1},'threshold'); assert(~exist('threshold','var'),'You cannot define globalData twice.'); threshold = varargin{2*i}; end;
    if ~ismember(varargin{2*i-1},[required, optional]); warning(['I could not identify the option ''' varargin{2*i-1} ''' and I am going to ignore it.']); end;
end

oldKnots = model.A.basis.knots;
filter = [boolean(ones(1, model.A.basis.degree + 1)), (slackVars' > threshold), boolean(ones(1, model.A.basis.degree + 1))]';
newKnots = oldKnots(filter);
T = [model.A model.B; model.C model.D];
basisNew = BSplineBasis(newKnots, T.basis.degree);
Tnew = project_to(T, basisNew);

Aproj = slice(Tnew,1:model.A.size(1),1:model.A.size(2));
Bproj = slice(Tnew,1:model.B.size(1),model.A.size(2)+1:model.A.size(2)+model.B.size(2));
Cproj = slice(Tnew,model.A.size(1)+1:model.A.size(1)+model.C.size(1),1:model.C.size(2));
Dproj = slice(Tnew,model.A.size(1)+1:model.A.size(1)+model.D.size(1),model.C.size(2)+1:model.C.size(2)+model.D.size(2));
regMdlProj = SSmod(Aproj,Bproj,Cproj,Dproj,schParam,model.Ts);

end

