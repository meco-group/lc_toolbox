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

function meas  = importBallscrewData(label, path, inputSignal, varargin)
%IMPORTBALLSCREWDATA Imports data in NetCDF format and casts it into a
%MeasurementData object.

    if ~isa(label,'char') || ~isa(path,'char') 
        error('You provided wrong input arguments: both the label and the path should be strings.'); 
    else if ~isa(inputSignal,'TimeSignal')
            error('Your input signal must be a TimeSignal object.');
         end
    end
    
    info = '';
    if nargin > 3
        if ~isa(varargin{1},'char'); error('Additional information should be provided in string format.'); end;
        info = varargin{1};
    else if nargin > 4
            error('Too many input arguments.');
         end
    end
    
    refTorque = ncread(path,'controller.ControlValues.0');
    refVMotor = ncread(path,'controller.ControlValues.2');
    refXMotor = ncread(path,'controller.ControlValues.1');
    measTorque = ncread(path,'plant.Measurements.2');
    measXLoad = ncread(path,'plant.Measurements.3');
    measXMotor = ncread(path,'plant.Measurements.1');
    measALoad = ncread(path,'plant.Measurements.4');
    
    timeStamp = ncread(path,'TimeStamp');
    
    try 
        file = dir(path);
        experimentDate = file.date;
    catch
        path = which(path);
        file = dir(path);
        experimentDate = file.date;
    end
    
    % we assume we work in torque mode
    meas = MeasurementData(label,inputSignal,[timeStamp refTorque measXMotor measXLoad measALoad],{'time','refTorque','measXMotor','measXLoad','measALoad'},experimentDate,info);

end

