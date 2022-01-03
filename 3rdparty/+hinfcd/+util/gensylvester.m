function [X,Y] = gensylvester(A,B,C,D,E,F,iterative)
% GENSYLVESTER Solves the generalized Sylvester equation
%
% [X,Y] = GENSYLVESTER(A,B,C,D,E,F) returns the solution to the generalized
% Sylvester equation
%   AX + YB = C
%   DX + YE = F
%
% [X,Y] = GENSYLVESTER(A,B,C,D,E,F,1) returns the same as the previous
% syntax, but solves the linear system iteratively instead of factorizing
% the system for backsubstitution

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
    assert(size(A,1)==size(A,2) && size(D,1)==size(D,2) && size(A,1)==size(D,1), 'GENSYLVESTER: A and D should be square matrices with the same size.'); 
    assert(size(B,1)==size(B,2) && size(E,1)==size(E,2) && size(B,1)==size(E,1), 'GENSYLVESTER: B and E should be square matrices with the same size.');
    assert(size(C,1)==size(A,1) && size(C,2)==size(B,2), 'GENSYLVESTER: C should have as many rows as A and as many columns as B.'); 
    assert(size(F,1)==size(D,1) && size(F,2)==size(E,2), 'GENSYLVESTER: F should have as many rows as F and as many columns as E.');
    m = size(A,1); n = size(B,1);
    
    % solve problem
    if nargin>6 && iterative
        % lsqr uses relative residual -> not always within the tolerance in the check below
        sol = lsqr([kron(eye(n),A) kron(B',eye(m)) ; kron(eye(n),D) kron(E',eye(m))], [C(:) ; F(:)], 1e-6, 200);
    else
        % the problem might be undetermined -> switch off warning
        warning('off','MATLAB:rankDeficientMatrix');
        sol = [kron(eye(n),A) kron(B',eye(m)) ; kron(eye(n),D) kron(E',eye(m))]\[C(:) ; F(:)];
        warning('off','MATLAB:rankDeficientMatrix');
    end

    % return solution
    X = reshape(sol(1:m*n),[m,n]);
    Y = reshape(sol(m*n+1:end),[m,n]);
    
    % check solution
    res = [A*X+Y*B-C ; D*X+Y*E-F];
    tol = 1e-6; 
    relres = norm(res)/norm([C(:) ; F(:)]);
    if relres>tol
        error(['GENSYLVESTER: Relative residual ' num2str(relres) ' is larger than the tolerance ' num2str(tol) '. Is the problem overdetermined or inconsistent?']); 
    end
    if any(isnan(res(:))) || any(isinf(res(:)))
        error('GENSYLVESTER: Relative residual contains NaN and/or Inf. Is the problem ill-conditioned or inconsistent?');
    end
    
end