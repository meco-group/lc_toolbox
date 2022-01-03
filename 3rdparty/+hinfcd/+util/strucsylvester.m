function [X,Y1] = strucsylvester(A,B,C,D,E,F,iterative)
% STRUCSYLVESTER Solves a specific generalized Sylvester equation
%
% [X,Y1] = STRUCSYLVESTER(A,B,C,D,E,F) returns the solution to the 
% generalized Sylvester equation with specific structure
%   A*[X1;X2] + [Y1;0]*B = C
%   D*X1 + Y1*E = F
% If D has as many rows as A, the command returns the same as
% GENSYLVESTER(A,B,C,D,E,F). If D, E and F have no rows, Y1 is empty and X
% is the solution of AX = C. 
%
% [X,Y1] = STRUCSYLVESTER(A,B,C,D,E,F,1) returns the same as the previous
% syntax, but solves the linear system iteratively instead of factorizing
% the system for backsubstitution
%
% See also HINFCD.UTIL.GENSYLVESTER.

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

    % check dimensions
    assert(size(D,1)==size(D,2), 'STRUCSYLVESTER: D should be a square matrix.'); 
    assert(size(B,1)==size(B,2) && size(E,1)==size(E,2), 'STRUCSYLVESTER: B and E should be square matrices.');
    assert(size(C,1)==size(A,1) && size(C,2)==size(B,2), 'STRUCSYLVESTER: C should have as many rows as A and as many columns as B.'); 
    assert(size(F,1)==size(D,1) && size(F,2)==size(E,2), 'STRUCSYLVESTER: F should have as many rows as D and as many columns as E.');
    n = size(B,1); nx = size(D,1); nu = size(A,2)-nx; nz = size(A,1)-nx;
    
    % solve problem
    UL = kron(eye(n),A);
    UR = kron(B', [eye(nx) ; zeros(nz,nx)]);
    BL = kron(eye(n),[D zeros(nx,nu)]);
    BR = kron(E',eye(nx)); 
    
    if nargin>6 && iterative
        % lsqr uses relative residual -> not always within the tolerance in the check below
        sol = lsqr([UL UR ; BL BR],[C(:) ; F(:)], 1e-6, 200);
    else
        % the problem might not have a unique solution -> switch off warning
        warning('off','MATLAB:rankDeficientMatrix');
        sol = [UL UR ; BL BR]\[C(:) ; F(:)];
        warning('on','MATLAB:rankDeficientMatrix');
    end
    
    % return solution
    X = reshape(sol(1:(n*(nx+nu))),[(nx+nu),n]);
    Y1 = reshape(sol((n*(nx+nu)+1):end),[nx,n]);
    
    % check solution
    res = [A*X+[Y1 ; zeros(nz,n)]*B-C ; D*X(1:nx,:)+Y1*E-F];
    tol = 1e-6; 
    relres = norm(res)/norm([C(:) ; F(:)]);
    if relres>tol
        error(['STRUCSYLVESTER: Relative residual ' num2str(relres) ' is larger than the tolerance ' num2str(tol) '. Is the problem overdetermined or inconsistent?']); 
    end
    if any(isnan(res(:))) || any(isinf(res(:)))
        error('STRUCSYLVESTER: Relative residual contains NaN and/or Inf. Is the problem ill-conditioned or inconsistent?');
    end
    
end