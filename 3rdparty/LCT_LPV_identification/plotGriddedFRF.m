function varargout = plotGriddedFRF( varargin )

required = {'griddedFRF'};
optional = {'oiChannel','figNumber','color'};

calls = varargin(1:2:end);
assert(mod(nargin,2) == 0,'Odd number of arguments.');
assert(iscellstr(calls),'The names of your name/value pairs should be strings.');
assert(~any(~ismember(required,calls)),['These arguments are required: ' strjoin(required)]);
for i = 1:nargin/2
    if strcmp(varargin{2*i-1},'griddedFRF'); assert(~exist('griddedFRF','var'),'You cannot define the griddedFRF twice.'); griddedFRF = varargin{2*i}; end;
    if strcmp(varargin{2*i-1},'oiChannel'); assert(~exist('oiChannel','var'),'You cannot define the oiChannel twice.'); oiChannel = varargin{2*i}; end;
    if strcmp(varargin{2*i-1},'figNumber'); assert(~exist('figNumber','var'),'You cannot define the figNumber twice.'); figNumber = varargin{2*i}; end;
    if strcmp(varargin{2*i-1},'color'); assert(~exist('color','var'),'You cannot define the color twice.'); color = varargin{2*i}; end;
    if ~ismember(varargin{2*i-1},[required, optional]); warning(['I could not identify the option ''' varargin{2*i-1} ''' and I am going to ignore it.']); end;
end
assert(isa(griddedFRF,'Gridmod'),'Unsupported data type.')

[mags, phss] = deal(zeros(numel(griddedFRF.grid_),numel(griddedFRF.grid_{1}.Frequency)));
[Xid,Yid] = meshgrid(griddedFRF.params_{2}', griddedFRF.grid_{1}.Frequency);      

if exist('oiChannel','var')
    assert(isnumeric(oiChannel),'Variable oiChannel must be a numerical array (vector).')  
else
    oiChannel = [1 1];
end
if exist('color','var')
    assert(ischar(color),'Variable color must be a character.')  
else
    color = 'b';
end

for k = 1:numel(griddedFRF.grid_)
    mags(k,:) = db(abs(squeeze(griddedFRF.grid_{k}.frdata(oiChannel(1),oiChannel(2),:))));
    phss(k,:) = (unwrap(phase(squeeze(griddedFRF.grid_{k}.frdata(oiChannel(1),oiChannel(2),:)))))*180/pi;
end

if exist('figNumber','var')
    assert(isnumeric(figNumber),'Variable figNumber must be a numerical array (vector).')  
    assert(numel(figNumber)==2,'Variable figNumber must have two elements.')  
    figure(figNumber(1))
    varargout{1} = plot3(Xid,Yid, mags',color,'linewidth', 1); grid on; hold on;
    ylabel('Frequency [Hz]');
    xlabel(griddedFRF.params_{1});
    zlabel('Magnitude [dB]');
    set(gca,'yscale','log');
    xlim([min(griddedFRF.params_{2}) max(griddedFRF.params_{2})])
    ylim([min(griddedFRF.grid_{1}.Frequency) max(griddedFRF.grid_{1}.Frequency)])
    view(125,50);
    figure(figNumber(2))
    varargout{2} = plot3(Xid, Yid, phss',color,'linewidth', 1); grid on; hold on;
    ylabel('Frequency [Hz]');
    xlabel(griddedFRF.params_{1});
    zlabel('Phase [^{o}]');
    set(gca,'yscale','log');
    xlim([min(griddedFRF.params_{2}) max(griddedFRF.params_{2})])
    ylim([min(griddedFRF.grid_{1}.Frequency) max(griddedFRF.grid_{1}.Frequency)])
    view(125,50);
else
    figure
    varargout{1} = plot3(Xid,Yid, mags',color,'linewidth', 1); grid on; hold on;
    ylabel('Frequency [Hz]');
    xlabel(griddedFRF.params_{1});
    zlabel('Magnitude [dB]');
    set(gca,'yscale','log');
    xlim([min(griddedFRF.params_{2}) max(griddedFRF.params_{2})])
    ylim([min(griddedFRF.grid_{1}.Frequency) max(griddedFRF.grid_{1}.Frequency)])
    view(125,50);
    figure
    varargout{2} = plot3(Xid, Yid, phss',color,'linewidth', 1); grid on; hold on;
    ylabel('Frequency [Hz]');
    xlabel(griddedFRF.params_{1});
    zlabel('Phase [^{o}]');
    set(gca,'yscale','log');
    xlim([min(griddedFRF.params_{2}) max(griddedFRF.params_{2})])
    ylim([min(griddedFRF.grid_{1}.Frequency) max(griddedFRF.grid_{1}.Frequency)])
    view(125,50);
end

end

