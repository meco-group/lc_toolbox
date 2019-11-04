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

function s = SUMmod(in1,in2)
%SUMmod Sum relation between signals
%   Defines the algebraic sum of 2, resulting in a 
%   new signal.

N = length(in1);
if length(in2) ~= N
    error('Inputs for SUM must be of equal length');
end
    
in = [in1;in2];
out = Signal(N);
A = zeros(0,0); B = zeros(0,2*N);
C = zeros(N,0); D = [eye(N),eye(N)];

s = SSmod(A,B,C,D,in,out);
end

