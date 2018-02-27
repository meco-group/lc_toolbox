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

function [IA,LOCB] = ismembertols(A,B,tol)
%ISMEMBERTOLS Summary of this function goes here
%   Detailed explanation goes here

if nargin < 3
    if eps < 1e-10
        tol = 1e-12; %double precision machine
    else
        tol = 1e-6;
    end
end

IA = false(1,length(A));
LOCB = zeros(1,length(A));

ib = 1;
for ia = 1:length(A)
    while B(ib) < (A(ia) - tol) 
        ib = ib+1;
        if ib>length(B), break; end
    end
    if ib>length(B), break; end
    if abs(B(ib)-A(ia)) < tol
        IA(ia) = true;
        LOCB(ia) = ib;
    end
end