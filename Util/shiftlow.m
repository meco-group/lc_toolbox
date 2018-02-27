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

function [sys_shift] = shiftlow(sys,varargin)
%SHIFTLOW Shifts low frequency poles to the origin
%   Detailed explanation goes here
if(nargin>1)
    pmin = varargin{1};
else
    pmin = 1e-2;
end

flag = false;
if ~isa(sys,'numlti')
    sys = std(sys);
    flag = true;
end

if ~all(size(sys)==1)
    for i = 1:size(sys,1)
        for k = 1:size(sys,2)
            sys(i,k) = shiftlow(sys(i,k),pmin);
        end
    end
    sys_shift = sys;
else
    [z,p,k] = zpkdata(sys,'v');

    p_shift = arrayfun(@(x) shift(x,pmin), p);
    z_shift = arrayfun(@(x) shift(x,pmin), z);
    
    [~,Ip] = sort(abs(p_shift),'ascend');
    p_shift = p_shift(Ip);
    [~,Iz] = sort(abs(z_shift),'ascend');
    z_shift = z_shift(Iz);
    
    while(z_shift(1) == 0)
        p_shift(1) = [];
        z_shift(1) = [];
    end
        
    sys_shift = zpk(z_shift,p_shift,k);
    sys_shift = sys_shift.*(abs(evalfr(sys,j*1e5))./abs(evalfr(sys_shift,j*1e5)));
    
%     if(angle(evalfr(sys_shift,j*1e5) - angle(evalfr(sys,j*1e5))
    
end

if flag
    sys_shift = fromstd(sys_shift);
end

end

function p = shift(p,varargin)
    if(nargin>1)
        pmin = varargin{1};
    else
        pmin = 1e-2;
    end
    p(abs(p)<pmin) = 0;
end