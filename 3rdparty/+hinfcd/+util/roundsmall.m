function A = roundsmall(A,zerotol)
% ROUNDSMALL Sets small values to zero
% Sets values of the array A that are smaller than a tolerance to zero
%
% Ar = ROUNDSMALL(A, tol) returns Ar, such that Ar(i) = A(i) if abs(A(i))>=
% zerotol, else Ar(i) = 0. The default value of zerotol is sqrt(eps). 

% This file is part of hinfcd.
% Copyright (c) 2019, Laurens Jacobs, MECO Research Team @ KU Leuven. 
% 
% hinfcd is free software: you can redistribute it and/or modify it under 
% the terms of the GNU Lesser General Public License as published by the 
% Free Software Foundation, version 3.
% 
% hinfcd is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for
% more details.
% 
% You should have received a copy of the GNU Lesser General Public License
% along with hinfcd. If not, see <https://www.gnu.org/licenses/>.

    if nargin<2
        zerotol = sqrt(eps);
    end
    A(abs(A)<zerotol) = 0;

end