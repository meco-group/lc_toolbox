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

classdef Order
    %ORDER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access=private)
        n = 0;
        s = [];
    end
    
    methods
        function self = Order(n,varargin)
            self.n = n;
            if nargin > 1
                self.s = varargin{1};
            end
        end
        
        function n = order(self)
            n = self.n;
        end
        
        function b = hassystem(self)
            b = ~isempty(self.s);
        end
        
        function s = system(self)
            s = self.s;
        end
    end
    
end

