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

classdef (InferiorClasses = {?zpk,?tf,?ss,?frd}) AnalyticModel < Model
    %ANALYTICMODEL Summary of this class goes here
    %   Detailed explanation goes here
    
    methods

        function self = AnalyticModel()
            self@Model();
        end
   
        function varargout = subsref(self,s)
            % SUBSREF reimplement subsref for systems subindexing            
            if strcmp(s(1).type,'()')
                siz = size(self);
                assert(~(isa(s(1).subs{1},'Signal') || isa(s(1).subs{2},'Signal')),'Cannot index a model using signals. Signals refer only to systems');
                if s(1).subs{1}==':';s(1).subs{1}=1:siz(1);end
                if s(1).subs{2}==':';s(1).subs{2}=1:siz(2);end
                assert(all(cellfun(@isnumeric,s(1).subs)),'only numeric index');
                assert(all(s(1).subs{1}<=siz(1)) && all(s(1).subs{2}<=siz(2)),'Index exceeds system dimensions.');
                varargout = {submodel(self,s(1).subs{1},s(1).subs{2})};
                
                if length(s)>1
                    varargout = {builtin('subsref',varargout{1},s(2:end))};
                end
            else
                varargout = {builtin('subsref',self,s(:))};
            end
        end
        
    end

end