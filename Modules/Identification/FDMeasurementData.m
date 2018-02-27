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

classdef FDMeasurementData < MeasurementData & FRDmod
    %MEASUREMENTDATA Class -- Summary of this class goes here
    %   Detailed explanation goes here
    
    methods
        
        function self = FDMeasurementData(varargin)
            
            idxlabel = [];
            idxdate = [];
            idxinfo = [];
            
            for i = 1:nargin/2
                if strcmp(varargin{2*i-1},'label'); assert(~exist('label','var'),'You cannot define a label twice.'); label = varargin{2*i}; idxlabel = i; end;
                if strcmp(varargin{2*i-1},'date'); assert(~exist('date','var'),'You cannot define a date twice.'); date = varargin{2*i}; idxdate = i; end;
                if strcmp(varargin{2*i-1},'info'); assert(~exist('info','var'),'You cannot define additional information twice.'); info = varargin{2*i}; idxinfo = i; end;
            end
            
            varargin([2*idxlabel,2*idxdate,2*idxinfo,2*idxlabel-1,2*idxdate-1,2*idxinfo-1]) = [];
            
            assert(exist('label','var')==1, 'Your measurement data should have a label.');
            if ~exist('date','var'); date = datestr(today); end;
            if ~exist('info','var'); info = ''; end;
            
            self@MeasurementData(label, info, date);
            self@FRDmod(varargin{:});
        end
        
        function display(self)
            if ndims(self) <= 3
                num = 1;
            else 
                num = size(self,4);
            end
            disp(['Frequency domain measurement data ''' self.label_ ''' recorded on ' self.date_ ' containing ' num2str(num) ' frequency response function(s) representing one excitation signal.']);
        end
        
        function disp(self)
            display(self)
        end
        
        function sys = std(self)
            p = properties(self);
            sys = frd(self.ResponseData,self.Frequency);
            for i = 1:length(p)
                if all(~strcmp(p{i},{'label_', 'info_', 'date_'}));
                    sys.(p{i}) = self.(p{i});
                end
            end
        end

    end

end
    

