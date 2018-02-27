% This file is part of LCToolbox.
% (c) Copyright 2018 - MECO Research Team, KU Leuven. 
%
% LCToolbox is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published 
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% LCToolbox is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU Lesser General Public License for more details.
% 
% You should have received a copy of the GNU Lesser General Public License
% along with LCToolbox. If not, see <http://www.gnu.org/licenses/>.

function [ out ] = importELMO( varargin )
% e.g. importELMO('label','ELMOexp_100%','date','12-June-2017','info','small mass added','dataPackName1','dataPackName2', ... , 'dataPackNameNdata', cutOpt, mergeOpt)
% cuttingOptions = [fmin1 fmax1; fmin2 fmax2; ...; fminNdata fmaxNdata]
% mergingOptions = {[dataPackIndex1 dataPackIndex2 ... ], [ ... ], [ ... ], ...}

infoCount = 0;
for i = 1:3;
    if strcmp(varargin{2*i-1},'label'); assert(~exist('label','var'),'You cannot define the label twice.'); label = varargin{2*i}; infoCount = infoCount+1; end;
    if strcmp(varargin{2*i-1},'date'); assert(~exist('date','var'),'You cannot define the date twice.'); date = varargin{2*i}; infoCount = infoCount+1; end;
    if strcmp(varargin{2*i-1},'info'); assert(~exist('info','var'),'You cannot define additional information twice.'); info = varargin{2*i}; infoCount = infoCount+1; end;
end
assert(exist('label','var')==1, 'Your measurement data should have a label.');
if ~exist('date','var'); date = datestr(today); end;
if ~exist('info','var'); info = ''; end;

optCount = 0;
for i = 1:2:3
    if strcmp(varargin{end-i},'cut'); assert(~exist('cut','var'),'You cannot define the cutting matrix twice.'); assert(isnumeric(varargin{end-i+1}),'Cutting matrix should be a numeric array.'); cutOpt = varargin{end-i+1}; optCount = optCount+1; end; 
    if strcmp(varargin{end-i},'merge'); assert(~exist('merge','var'),'You cannot define the merging matrix twice.'); assert(iscell(varargin{end-i+1}),'Merge options should be a cell array.'); mergeOpt = varargin{end-i+1}; optCount = optCount+1; end; 
end
ndata = nargin-2*infoCount-2*optCount;
 
% preallocation
load(varargin{2*infoCount+ndata})
dataStruct(ndata).Frequency = VelPosPlants(numel(VelPosPlants)).Freqs;  
dataStruct(ndata).Response = db2mag(VelPosPlants(numel(VelPosPlants)).Magnitude).*(cos(VelPosPlants(numel(VelPosPlants)).Phase*pi/180) + 1i*sin(VelPosPlants(numel(VelPosPlants)).Phase*pi/180));  
Ts = GantryController.SampleTime*1e-6;

% loading
for i = (2*infoCount+1):(2*infoCount+ndata-1)
    assert(ischar(varargin{i}),'The data file name should be a string.');
    try
        load(varargin{i})
    catch
        error('File %s not found.', varargin{i})
    end
    dataStruct(i-2*infoCount).Frequency = VelPosPlants(numel(VelPosPlants)).Freqs;
    dataStruct(i-2*infoCount).Response = db2mag(VelPosPlants(numel(VelPosPlants)).Magnitude).*(cos(VelPosPlants(numel(VelPosPlants)).Phase*pi/180) + 1i*sin(VelPosPlants(numel(VelPosPlants)).Phase*pi/180)); 
end

% cutting
if exist('cutOpt','var') 
    assert(size(cutOpt,1) == ndata, 'If provided, cutting options should have the number of rows equal to the number of data packages.')
    assert(size(cutOpt,2) == 2, 'If provided, cutting options should have two columns, one for the minimal frequencies, one for the and maximal.');
    for i = 1:ndata
        [~,indexFmin] = min(abs(dataStruct(i).Frequency - cutOpt(i,1)));
        [~,indexFmax] = min(abs(dataStruct(i).Frequency - cutOpt(i,2)));
        dataStruct(i).Frequency = dataStruct(i).Frequency(indexFmin:indexFmax);  
        dataStruct(i).Response = dataStruct(i).Response(indexFmin:indexFmax);  
    end
else
    for i = 1:ndata
        stop = find(dataStruct(i).Frequency == 0);
        dataStruct(i).Frequency = dataStruct(i).Frequency(1:(stop-1));  
        dataStruct(i).Response = dataStruct(i).Response(1:(stop-1));
    end
end   

% merging
if exist('mergeOpt','var')
    if ~isempty(mergeOpt)
        dataStructNew(numel(mergeOpt)).Frequency = [];
        dataStructNew(numel(mergeOpt)).Response = [];
        for i = 1:numel(mergeOpt)
            dataAcc = dataStruct(mergeOpt{i}(1));
            if numel(mergeOpt{i}) < 2
                warning('At least two elements needed for a merge.')
            else
                for j = 2:numel(mergeOpt{i}) 
                    dataAcc = merge2frfs(dataAcc, dataStruct(mergeOpt{i}(j)));
                end
            end
            dataStructNew(i) = dataAcc;
        end
    else
        dataStructNew = dataStruct;
    end
else
    dataStructNew = dataStruct;
end

% output
out = cell(1,numel(dataStructNew));
for i = 1:numel(dataStructNew)
    out{i} = FDMeasurementData('label',label,'date',date,'info',info,dataStructNew(i).Response, dataStructNew(i).Frequency, Ts, 'FrequencyUnit', 'Hz', 'InputName', 'current', 'OutputName', 'velocity');
end

end

