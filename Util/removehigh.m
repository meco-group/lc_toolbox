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

function [sys_remove] = removehigh(sys,varargin)
%REMOVEHIGH Summary of this function goes here
%   Detailed explanation goes here
if(nargin>1)
    if(length(varargin{1}) == 1)
        xmin = [varargin{1};+inf];
    else
        xmin = varargin{1};
    end
else
    xmin = [1e5;inf];
end

flag = false;
if ~isa(sys,'numlti')
    stdsys = std(sys);
    flag = true;
else
    stdsys = sys;
end

[z,p,k] = zpkdata(stdsys);

p_rm = cellfun(@(x) remove(x,xmin), p, 'UniformOutput', false);
z_rm = cellfun(@(x) remove(x,xmin), z, 'UniformOutput', false);

sys_remove = zpk(z_rm,p_rm,k);
sys_remove = sys_remove.*(abs(evalfr(stdsys,j*1e-2))./abs(evalfr(sys_remove,j*1e-2)));
sys_remove = correctPhase(sys_remove,stdsys);

if flag
    sys_remove = fromstd(sys_remove);
    sys_remove.name = sys.name;
end

end

function x = remove(x,varargin)
    if(nargin>1)
        if(length(varargin{1}) == 1)
            xmin = [varargin{1};+inf];
        else
            xmin = varargin{1};
        end
    else
        xmin = [1e5;+inf];
    end
    x((abs(x)>xmin(1)) & (abs(x)<xmin(2))) = [];
end

function [sys_remove] = correctPhase(sys_remove,sys)
    phase = angle(evalfr(sys_remove,0));
    phaseRef = angle(evalfr(sys,0));
    correction = arrayfun(@(x,y) (abs(x-y) < pi/2), angle(evalfr(sys_remove,j*1e-2)), angle(evalfr(sys,j*1e-2)));
    correction = double(correction);
    correction(correction == 0) = -1;
    
    sys_remove = sys_remove.*correction;
end