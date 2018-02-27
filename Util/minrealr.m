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

function [sys] = minrealr(sys,varargin)
%MINREALR Relative minreal
%   Cancel out poles and zeros based on relative distance rather than
%   absolute distance. This is currently implemented for siso only

if(~all(size(sys)==1))
%     for k = 1:size(sys,1)
%         for l = 1:size(sys,2)
%             sys(k,l) = minrealr(sys(k,l),varargin{:});
%         end
%     end
else
    [z,p,k] = zpkdata(sys,'v');
    tol = 1e-3;
    if(nargin>1)
        tol = varargin{1};
    end

    remove_z = false(size(z));
    remove_p = false(size(p));

    for(j=length(z):-1:1)
        for(i = length(p):-1:1)
            az = abs(z(j)); ap = abs(p(i));
            if(2*(abs(az-ap))/(az+ap) < tol)
%                 disp('removed pole/zero');
%                 disp('pole'); disp(p(i));
%                 disp('zero'); disp(z(j));  
                z(j) = []; p(i) = [];   
                break;
            end
        end
    end

    sys = zpk(z,p,k,sys.Ts);
end
end

