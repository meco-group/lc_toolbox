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

classdef SchedulingParameter < splines.Function
    properties
        rate_
        bounded
    end
    
    methods
        function parameter = SchedulingParameter(argument, range, rate)
            if isnumeric(range)
                b = splines.BSplineBasis(range,1,2);
                c = range';
            else
                b = range.basis;
                c = range.coeff;
            end
            parameter@splines.Function(splines.TensorBasis(b, argument), c);
            parameter.bounded = 0;            
            switch nargin
                case 3
                    if ~all(rate==0)
                        parameter.bounded = 1;
                    else
                        parameter.bounded = 0;
                    end

                    if length(rate) == 1
                        rate = [-abs(rate), abs(rate)];
                    elseif length(rate) == 2
                        rate = sort(rate);
                    else
                        error('rate should have one or two dimensions');
                    end
                    parameter.rate_ = rate;
                case 2
                    parameter.rate_ = [0, 0];
                otherwise
                    error(['SchedulingParameter expect 2 or 3 inputs: argument(string), range(interval), [rate(interval)], not compatible with ', num2str(nargin) ,'provided']);
            end
        end
        
        function s = to_string(self)
            s = [to_string@splines.Function(self), char(10), 'and a rate of variation bounded by ', num2str(self.rate_)];
        end
        
        function r = range(self)
            r = [self.basis.domain.min, self.basis.domain.max];
        end
        
        function r = rate(self)
            r = self.rate_;
        end
        
        function r = add_parameter_derivative(self,c_or_d)
            r = {};
             assert(self.bounded==1,'cannot add, parameter is not bounded')
                    switch c_or_d
                        case 'c'
                            p = SchedulingParameter(strcat(self.tensor_basis.arguments,'dot'),self.rate,0);
                            r = {self,p};
                        case 'd'
                            p = SchedulingParameter(strcat('d',self.tensor_basis.arguments),self.rate,0);
                            r = {self,p};
                    end
        end

        function a = argument(self)
            a = arguments(tensor_basis(self));
        end
end    
    methods(Static)
        function [LIA,LOCB] = ismember(A,B)
            argsA = cellfun(@(x)arguments(tensor_basis(x)),A,'un',0);
            argsB = cellfun(@(x)arguments(tensor_basis(x)),B,'un',0);
            [LIA,LOCB] = ismember(argsA,argsB);
        end
    end
    

end