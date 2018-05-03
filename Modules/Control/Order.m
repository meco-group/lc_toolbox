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

classdef Order < Specification
    % Order is a specific Specification type defining the desired
    % controller(s) order(s).
    
    properties (Access=private)
        n = 0;      % desired order @type double
        s = [];     % system to which the order constraint applies @type AbstractSystem
    end
    
    methods
        function self = Order(n,varargin)
        % Constructor for Order objects.
        % 
        % Parameters:
        %   n: the desired order @type double
        %   varargin: can contain a system in case the order constraint only holds for a specific (sub)system that is optimized for @type AbstractSystem
        %
        % Return values:
        %   self: the order specification @type Order
        
            assert(mod(n,1)==0, 'Orders should have an integer number.');
            self.n = n;
            if nargin > 1
                self.s = varargin{1};
            end
        end
        
        function n = order(self)
        % Returns the specified order. 
        %
        % Return values:
        %   n: the specified order @type double
            n = self.n;
        end
        
        function b = hassystem(self)
        % Checks whether the order constraint applies to a specific
        % (sub)system.
        %
        % Return values:
        %   b: true if the order constraint applies to a specific
        %   (sub)system @type logical
            b = ~isempty(self.s);
        end
        
        function s = system(self)
        % Returns the (sub)system to which the order constraint applies.
        %
        % Return values:
        %   s: the (sub)system to which the order constraint applies @type
        %   AbstractSystem
            s = self.s;
        end
    end
    
end

