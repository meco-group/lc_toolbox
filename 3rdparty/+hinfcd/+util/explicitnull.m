function A = explicitnull(E,nofcol)
% EXPLICITNULL Explicit formulation of equality constraints
%
% A = EXPLICITNULL(E,nofcol) reparametrizes the equality constraint 
% E'*X = 0 as a constraint vec(X) = A(:,1)*x1 + A(:,2)*x2 + ... where the 
% optimization matrix variable X is transformed to an optimization vector 
% variable [x1 x2 ...]'. The length of this vector depends on the rank of
% E. E is assumed to be square, X is not necessarily. X has nofcol columns.
%
% See also HINFCD.UTIL.EXPLICITSYM. 
%
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

   B = kron(eye(nofcol),E'); 
   A = null(B); 
   
end