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

classdef (InferiorClasses = {?zpk,?tf,?ss,?frd}) LTILFTmod < AbstractLFTmod & AbstractLTImod & AnalyticModel
    %LTILFTMOD Summary of this class goes here
    %   Detailed explanation goes here
    
    methods
        function self = LTILFTmod(M,Nu,Nl,E,Ts)
            if iscell(M), checkM = all(cellfun(@isnumeric,M(:)));
            else checkM = isnumeric(M); end
            assert(checkM&&isnumeric(Nu)&&isnumeric(Nl)&&isnumeric(E),'All inputs should be numeric');
            self = self@AbstractLFTmod(M,Nu,Nl,E,Ts);
        end
        
        function sys = std(self)
            S = transpose(lft2ss(self));
            sys = dss(S{:},self.E,self.Ts);
        end
        
        function sys = simplify(self)
            S = transpose(lft2ss(self));
            sys = LTIDSSmod(S{:},self.E,self.Ts);
        end
        
        function product = mtimes(self,other)
            if isnumeric(self) || isnumeric(other)
                product = mtimes@AbstractLFTmod(self,other);
            else
                product = mtimes@AnalyticModel(self,other);
            end
        end
    end
end

