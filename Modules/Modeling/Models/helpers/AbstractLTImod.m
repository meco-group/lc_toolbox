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

classdef AbstractLTImod < AbstractLPVmod
% AbstractLTImod provides a blueprint consisting of the properties and
% methods every LTI model of LCToolbox should have. 
    
        
    methods (Abstract)
        std(self);
        % Transforms the LTI model into a standard MATLAB model.
    end
    
    methods        
        function self = AbstractLTImod
        % Constructor for AbstractLTImod objects. 
            self@AbstractLPVmod({});
        end
        
        %% tools
        
        function varargout = pole(varargin)
        % Returns the poles of the LTI model.
        % Overloads MATLAB's \c pole() for AbstractLTImod objects.
            varargin = stdargs(varargin);
            [varargout{1:nargout}] = pole(varargin{:});
        end
        
        function varargout = zero(varargin)
        % Returns the zeros of the LTI model.
        % Overloads MATLAB's \c zero() for AbstractLTImod objects.
            varargin = stdargs(varargin);
            [varargout{1:nargout}] = zero(varargin{:});
        end
                
        function b = isstable(varargin)
        % Checks the stability of the LTI model.
        % Overloads MATLAB's \c isstable() for AbstractLTImod objects.
            varargin = stdargs(varargin);
            b = isstable(varargin{:});
        end
        
        function n = norm(varargin)
        % Returns the norm of the LTI model.
        % Overloads MATLAB's \c norm() for AbstractLTImod objects.
            varargin = stdargs(varargin);
            n = norm(varargin{:});
        end
        
        function [GS,GNS] = stabsep(self)
        % Decomposes the LTI models into a stable and an unstable part. 
        % Overloads MATLAB's \c stapsep() for AbstractLTImod objects.
            if isempty(pole(self))
                GS = self;
                GNS = fromstd(zeros(size(self)));
            else
                [GS,GNS] = stabsep(std(self));
                GS = fromstd(GS); GNS = fromstd(GNS);
            end
        end
        
        function i = inv(self)
        % Inverts the LTI model.
        % Overloads MATLAB's \c stapsep() for AbstractLTImod objects.
            i = fromstd(inv(std(self)));
        end
        
        function varargout = lsim(varargin)
        % Simulates the LTI model in the time domain.
        % Overloads MATLAB's \c lsim() for AbstractLTImod objects.
            varargin = stdargs(varargin);
            [varargout{1:nargout}] = lsim(varargin{:});
        end
        
        function varargout = c2d(self,varargin)
        % Discretizes the continuous-time LTI model.
        % Overloads MATLAB's \c c2d() for AbstractLTImod objects.
            [sys,G] = c2d(std(self),varargin{:});
            varargout{1} = fromstd(sys);
            if nargout >= 2
                varargout{2} = G;
            end
        end
        
        function varargout = d2c(self,varargin)
        % Interpolates the discrete-time LTI model.
        % Overloads MATLAB's \c d2c() for AbstractLTImod objects.
            [sys,G] = d2c(std(self),varargin{:});
            varargout{1} = fromstd(sys);
            if nargout >= 2
                varargout{2} = G;
            end
        end
        
    end
end

