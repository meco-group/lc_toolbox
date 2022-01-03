function [A,x0] = explicitsym(E,X0)
% EXPLICITSYM Explicit formulation of symmetry constraints 
%
% A = EXPLICITSYM(E) reparametrizes the equality constraint E'*X = X'*E as 
% a constraint vec(X) = A(:,1)*x1 + A(:,2)*x2 + ... where the optimization
% matrix variable X is transformed to an optimization vector variable 
% [x1 x2 ...]'. The length of this vector depends on the rank of E. E and X
% are assumed to have the same size.
%
% [A,x0] = EXPLICITSYM(E,X0) reparametrizes the equality constraint 
% E'*[X0 X] = [X0 X]'*E, where X0 is a known constant part, as a constraint
% vec(X) = x0 + A(:,1)*x1 + A(:,2)*x2 + ... where the optimization
% matrix variable X is transformed to an optimization vector variable 
% [x1 x2 ...]'. The length of this vector depends on the rank of E. E and X
% are assumed to have the same size.
%
% See also HINFCD.UTIL.EXPLICITNULL. 
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

   [n,m] = size(E); 
   
   % elaborate constraints as a homogeneous linear system
   B = zeros(0.5*m*(m-1),n*m);
   nr = 0;
   for i=1:m
       B(nr+(1:m-i),(i-1)*n+(1:n)) = -E(:,i+1:m)';
       B(nr+(1:m-i),(i*n+1):end) = kron(eye(m-i),E(:,i)');
       nr = nr+(m-i);
   end
   
   % add constant part to the equations
   if nargin>1
       assert(size(X0,1)==n && size(X0,2)<=m, 'EXPLICITSYM: X0 cannot have more columns than E, and should have as many rows as E.'); 
       B = [B ; eye(n*size(X0,2)) zeros(n*size(X0,2),n*(m-size(X0,2)))]; 
       b = [zeros(0.5*m*(m-1),1) ; X0(:)];
       warning('off','MATLAB:rankDeficientMatrix');
       x0 = B\b; 
       warning('on','MATLAB:rankDeficientMatrix');
       tol = 1e-6;
       relres = norm(B*x0-b)/norm(b); 
       assert(relres<tol, 'EXPLICITSYM: Overdetermined system. Does X0 satisfy the symmetry constraint?'); 
   else
       x0 = zeros(n*m,1); 
   end
   
   % parametrize the solution using a basis of the kernel
   A = null(B);
  
end