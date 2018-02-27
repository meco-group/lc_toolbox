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

classdef (Abstract) MeasurementData 
    %MEASUREMENTDATA Class -- Summary of this class goes here
    %   Detailed explanation goes here
    
    properties %(Access = private)
        date_ = '';      % string
        label_ = '';     % label for the data set
        info_ = '';      % string 
    end
    
    methods
        function self = MeasurementData(label, info, date)
            self.label_ = label;
            self.info_ = info;
            self.date_ = date;
        end
    end
    
end