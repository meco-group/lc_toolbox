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

classdef (InferiorClasses = {?zpk,?tf,?ss,?frd}) LTIDSSmod < AbstractDSSmod & AbstractLTImod & Model
    %LTIDSSMOD Summary of this class goes here
    %   Detailed explanation goes here
    
    methods
        function self = LTIDSSmod(A,B,C,D,E,varargin)
            assert(all([isnumeric(A),isnumeric(B),isnumeric(C),isnumeric(D),isnumeric(E)]),'All inputs should be numeric');
            self = self@AbstractDSSmod(A,B,C,D,E,varargin{:});
        end
        
        function sys = std(self)
            sys = dss(self.A,self.B,self.C,self.D,self.E,self.Ts);
        end
        
        function sys = simplify(self)
            sys = self;
        end
        
        function product = mtimes(self,other)
            if isnumeric(self) || isnumeric(other)
                product = mtimes@AbstractDSSmod(self,other);
            else
                product = mtimes@Model(self,other);
            end
        end
        
        function varargout = c2d(self,varargin)
            % TO DO: make difference between builtin and own implementation
            [sys,G] = c2d(std(self),varargin{:});
            varargout{1} = fromstd(sys);
            if nargout >= 2
                varargout{2} = G;
            end
        end
        
        function varargout = d2c(self,varargin)
            [sys,G] = d2c(std(self),varargin{:});
            varargout{1} = fromstd(sys);
            if nargout >= 2
                varargout{2} = G;
            end
        end
        
    end
end