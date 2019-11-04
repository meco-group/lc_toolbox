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

function [sys_remove] = removelow(sys,varargin)
%REMOVELOW Summary of this function goes here
%   Detailed explanation goes here
if(nargin>1)
    xmax = varargin{1};
else
    xmax = 1e-2;
end

[z,p,k] = zpkdata(sys);

p_rm = cellfun(@(x) remove(x,xmax), p, 'UniformOutput', false);
z_rm = cellfun(@(x) remove(x,xmax), z, 'UniformOutput', false);

sys_remove = zpk(z_rm,p_rm,k);
sys_remove = sys_remove.*(abs(evalfr(sys,j*1e5))./abs(evalfr(sys_remove,j*1e5)));
sys_remove = correctPhase(sys_remove,sys);

end

function x = remove(x,varargin)
    if(nargin>1)
        xmax = varargin{1};
    else
        xmax = 1e-2;
    end
    x(abs(x)<xmax) = [];
end

function [sys_remove] = correctPhase(sys_remove,sys)
    phase = angle(evalfr(sys_remove,j*1e5));
    phaseRef = angle(evalfr(sys,j*1e5));
    correction = arrayfun(@(x,y) (abs(x-y) < pi/2), angle(evalfr(sys_remove,j*1e5)), angle(evalfr(sys,j*1e5)));
    correction = double(correction);
    correction(correction == 0) = -1;
    
    sys_remove = sys_remove.*correction;
end

