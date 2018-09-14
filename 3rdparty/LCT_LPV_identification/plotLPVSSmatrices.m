function plotLPVSSmatrices( varargin )

required = {'LPVmod'};
optional = {'figNumber','color'};

calls = varargin(1:2:end);
assert(mod(nargin,2) == 0,'Odd number of arguments.');
assert(iscellstr(calls),'The names of your name/value pairs should be strings.');
assert(~any(~ismember(required,calls)),['These arguments are required: ' strjoin(required)]);
for i = 1:nargin/2
    if strcmp(varargin{2*i-1},'LPVmod'); assert(~exist('LPVmod','var'),'You cannot define the LPVmod twice.'); LPVmod = varargin{2*i}; end;
    if strcmp(varargin{2*i-1},'figNumber'); assert(~exist('figNumber','var'),'You cannot define the figNumber twice.'); figNumber = varargin{2*i}; end;
    if strcmp(varargin{2*i-1},'color'); assert(~exist('color','var'),'You cannot define the color twice.'); color = varargin{2*i}; end;
    if ~ismember(varargin{2*i-1},[required, optional]); warning(['I could not identify the option ''' varargin{2*i-1} ''' and I am going to ignore it.']); end;
end
assert(isa(LPVmod,'LPVDSSmod'),'Unsupported data type.')

if exist('color','var')
    assert(ischar(color),'Variable color must be a character.')  
else
    color = 'b';
end

schedGrid = linspace(LPVmod.getdssdata.domain.domain.min,LPVmod.getdssdata.domain.domain.max,1e3);
T = [LPVmod.A LPVmod.B; LPVmod.C LPVmod.D];
T_eval = T.list_eval(schedGrid);
T_eval_knots = T.list_eval(T.basis.knots);

if exist('figNumber','var')
    assert(isnumeric(figNumber),'Variable figNumber must be have a numeric value.')    
    figure(figNumber)
    for i = 1:T.size(1)
        for j = 1:T.size(2)
            subplot(T.size(1),T.size(2),(i-1)*T.size(2)+j)
            plot(schedGrid,T_eval(:,i,j),color,'Linewidth',1.5); hold on
            scatter(T.basis.knots,T_eval_knots(:,i,j),'MarkerEdgeColor',color); hold on
        end
    end  
else
    figure
    for i = 1:T.size(1)
        for j = 1:T.size(2)
            subplot(T.size(1),T.size(2),(i-1)*T.size(2)+j)
            plot(schedGrid,T_eval(:,i,j),color,'Linewidth',1.5); hold on
            scatter(T.basis.knots,T_eval_knots(:,i,j),'MarkerEdgeColor',color); hold on
        end
    end  
end

end

