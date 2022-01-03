function B = orthtol(A,zerotol)
% ORTHTOL Returns an orthonormal basis of the range of a matrix
%
% Same functionality as b = orth(a) from MATLAB, but allows for a
% tolerance. 
%
% b = ORTHTOL(a) returns a matrix b such that its columns form an
% orthonormal basis of Col(a)
%
% b = ORTHTOL(a,zerotol) returns the same as the previous syntax, but
% treats the singular values that are smaller than zerotol (default value:
% sqrt(eps)) as exact zeros

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

    if nargin>0; zerotol = sqrt(eps); end
    [U,S] = svd(A,'econ'); 
    B = U(:,abs(diag(S))>=zerotol);

end