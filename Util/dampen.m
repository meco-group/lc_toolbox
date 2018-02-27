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

function [sys_damped] = dampen(sys, varargin)
    range = [0,+inf];
    damping = 1;
    
    if(nargin>1)
        if(length(varargin{1})>1)
            range = varargin{1};
            if(nargin>2)
                damping = varargin{2};
            end
        else
            damping = varargin{1};
            if(nargin>2)
                range = varargin{2};
            end
        end
    end
    
    [z,p,k] = zpkdata(sys,'v');
    
    iz = (abs(z)>range(1)) & (abs(z)<range(2));
    ip = (abs(p)>range(1)) & (abs(p)<range(2));
    
    z(iz) = abs(z(iz)).*(sign(real(z(iz)))*damping + j*sign(imag(z(iz)))*sqrt(1-damping^2));
    p(ip) = abs(p(ip)).*(sign(real(p(ip)))*damping + j*sign(imag(p(ip)))*sqrt(1-damping^2));
    sys_damped = zpk(z,p,k);
end