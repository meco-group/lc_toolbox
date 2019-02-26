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

classdef Region < Specification
    
    % This class will allow interactive specification of regions for pole
    % placement in the near future.
    properties
        L = [];
        M = [];
    end
    
    methods
        function self = Region(type,varargin)
            if nargin > 1
                param = varargin{1};
            end
            switch type
                case 'half'
                    assert(isreal(param.alphal),'you should enter the left most real desired pole');
                    if(~isfield(param,'alphar'))
                        param.alphar = 0;
                    end
                    self.L = [2*param.alphal 0; 0 -2*param.alphar];
                    self.M = [-1 0; 0 1];
                case 'conic'
                    self.L = [0 0; 0 0];
                    self.M = [cos(param.beta) -sin(param.beta); sin(param.beta) cos(param.beta)];
                otherwise
                    error('For now, it can only take the left half section plane or a conic plane including damping, denoted by either "half" or "conic"');
            end
        end
        function m = region(self)
        m.L = self.L;
        m.M = self.M;
        end
    end
end

