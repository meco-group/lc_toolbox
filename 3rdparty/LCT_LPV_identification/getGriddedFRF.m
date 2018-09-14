function varargout = getGriddedFRF( varargin )

required = {'LPVmod','frequencyGrid','frequencyUnit','schedParam','schedGrid'};
optional = {'measGriddedFRF'};

calls = varargin(1:2:end);
assert(mod(nargin,2) == 0,'Odd number of arguments.');
assert(iscellstr(calls),'The names of your name/value pairs should be strings.');
assert(~any(~ismember(required,calls)),['These arguments are required: ' strjoin(required)]);
for i = 1:nargin/2
    if strcmp(varargin{2*i-1},'LPVmod'); assert(~exist('LPVmod','var'),'You cannot define the LPVmod twice.'); LPVmod = varargin{2*i}; end;
    if strcmp(varargin{2*i-1},'frequencyGrid'); assert(~exist('frequencyGrid','var'),'You cannot define the frequencyGrid twice.'); frequencyGrid = varargin{2*i}; end;
    if strcmp(varargin{2*i-1},'frequencyUnit'); assert(~exist('frequencyUnit','var'),'You cannot define the color twice.'); frequencyUnit = varargin{2*i}; end;
    if strcmp(varargin{2*i-1},'schedParam'); assert(~exist('schedParam','var'),'You cannot define the schedParam twice.'); schedParam = varargin{2*i}; end;
    if strcmp(varargin{2*i-1},'schedGrid'); assert(~exist('schedGrid','var'),'You cannot define the schedGrid twice.'); schedGrid = varargin{2*i}; end;
    if strcmp(varargin{2*i-1},'measGriddedFRF'); assert(~exist('measGriddedFRF','var'),'You cannot define the measGriddedFRF twice.'); measGriddedFRF = varargin{2*i}; end;
    if ~ismember(varargin{2*i-1},[required, optional]); warning(['I could not identify the option ''' varargin{2*i-1} ''' and I am going to ignore it.']); end;
end
assert(isa(LPVmod,'LPVDSSmod'),'Unsupported data type.')
assert(isa(measGriddedFRF,'Gridmod'),'Unsupported data type.')
assert(ischar(frequencyUnit),'Variable frequencyUnit must be a string.')
assert(isa(schedParam,'SchedulingParameter'),'Variable schedParam must be of the SchedulingParameter class.')
assert(iscell(schedGrid)&&(numel(schedGrid)==2),'Variable schedGrid must be a cell array consisting of two elements.')

if exist('measGriddedFRF','var')
    if strcmp(measGriddedFRF.grid_{1}.FrequencyUnit,frequencyUnit)
        assert(isequal(measGriddedFRF.grid_{1}.Frequency,frequencyGrid),'The provided frequencyGrid and the frequency grid of the measured FRFs should coincide.')
    elseif strcmp(measGriddedFRF.grid_{1}.FrequencyUnit,'Hz') && strcmp(frequencyUnit,'rad/s')
        assert(isequal(2*pi*measGriddedFRF.grid_{1}.Frequency,frequencyGrid),'The provided frequencyGrid and the frequency grid of the measured FRFs should coincide.')
    elseif strcmp(measGriddedFRF.grid_{1}.FrequencyUnit,'rad/s') && strcmp(frequencyUnit,'Hz')
        assert(isequal(measGriddedFRF.grid_{1}.Frequency,2*pi*frequencyGrid),'The provided frequencyGrid and the frequency grid of the measured FRFs should coincide.')
    end
end

[mdlCells, mdlFRFcells, mdlErrFRFcells] = deal(cell(1,numel(schedGrid(2))));
if strcmp(frequencyUnit,'Hz')
    for i = 1:numel(schedGrid{2})
        mdlCells{i} = LPVmod.evalme(schedGrid{2}(i),{schedParam});
        mdlFRFcells{i} = FDMeasurementData('label', 'label1', 'date', '27-October-2017', freqresp(mdlCells{i},2*pi*frequencyGrid), frequencyGrid, LPVmod.Ts, 'FrequencyUnit', frequencyUnit, 'InputName', 'input1', 'OutputName', 'output1');
        if exist('measGriddedFRF','var')
            mdlErrFRFcells{i} = FDMeasurementData('label', 'label1', 'date', '27-October-2017', abs(measGriddedFRF.grid_{i}.frdata-freqresp(mdlCells{i},2*pi*frequencyGrid)), frequencyGrid, LPVmod.Ts, 'FrequencyUnit', frequencyUnit, 'InputName', 'input1', 'OutputName', 'output1');
        end
    end
elseif strcmp(frequencyUnit,'rad/s')
    for i = 1:numel(schedGrid{2})
        mdlCells{i} = LPVmod.evalme(schedGrid{2}(i),{schedParam});
        mdlFRFcells{i} = FDMeasurementData('label', 'label1', 'date', '27-October-2017', freqresp(mdlCells{i},frequencyGrid), frequencyGrid, LPVmod.Ts, 'FrequencyUnit', frequencyUnit, 'InputName', 'input1', 'OutputName', 'output1');
        if exist('measGriddedFRF','var')
            mdlErrFRFcells{i} = FDMeasurementData('label', 'label1', 'date', '27-October-2017', abs(measGriddedFRF.grid_{i}.frdata-freqresp(mdlCells{i},frequencyGrid)), frequencyGrid, LPVmod.Ts, 'FrequencyUnit', frequencyUnit, 'InputName', 'input1', 'OutputName', 'output1');
        end
    end
end

varargout{1} = Gridmod(mdlFRFcells,schedGrid);
if exist('measGriddedFRF','var')
    varargout{2} = Gridmod(mdlErrFRFcells,schedGrid);
end

end

