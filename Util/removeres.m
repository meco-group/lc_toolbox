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

function [sys_remove] = removeres(sys,w,tol)
%REMOVERES Removes the resonanses from a given system
%   Detailed explanation goes here

[z,p,k] = zpkdata(sys);

p_rm = cellfun(@(x) remove(x,w,tol), p, 'UniformOutput', false);
z_rm = cellfun(@(x) remove(x,w,tol), z, 'UniformOutput', false);

sys_remove = zpk(z_rm,p_rm,k);
sys_remove = sys_remove.*(abs(evalfr(sys,0)./abs(evalfr(sys_remove,0))));
sys_remove = correctPhase(sys_remove,sys);

end

function x = remove(x,w,tol)
    x((abs(real(x))<tol) & (abs(x) > w)) = [];
end

function [sys_remove] = correctPhase(sys_remove,sys)
    phase = angle(evalfr(sys_remove,0));
    phaseRef = angle(evalfr(sys,0));
    correction = arrayfun(@(x,y) (abs(x-y) < pi/2), angle(evalfr(sys_remove,0)), angle(evalfr(sys,0)));
    correction = double(correction);
    correction(correction == 0) = -1;
    
    sys_remove = sys_remove.*correction;
end

