function [obj,u,v,feigv] = kronreal(obj,zerotol)
% KRONREAL Kronecker-Weierstrass canonical realization of a descriptor system
% The Kronecker-Weierstrass form is found in two steps:
%   1. A triangular form, i.e. such that the pencil A - sE is upper
%   triangular. 
%   2. A block-diagonal form, i.e. such that the slow and the fast
%   subsystem are completely decoupled. 
% This procedure is based on M. Gerdin, "Computation of a canonical form
% for linear differential-algebraic equations", Technical report: 
% LiTH-ISY-R-2602, Linköping University, 2004.
%
% Note: The Kronecker-Weierstrass decomposition is a cumbersome realization
% from a numerical point of view, i.e., it is usually ill-conditioned. It
% is strongly encouraged to check the accuracy of the results thorougly
% and, if required, change the zerotol value (see below). 
%
% b = KRONREAL(a) returns the Kronecker-Weierstrass realization b of a
%
% [b,u,v] = KRONREAL(a) returns the same as the previous syntax, and two
% matrices u and v such that (u*E*v, u*A*v, u*B, C*v) is the
% Kronecker-Weierstrass decomposition of (E,A,B,C) of a
%
% [b,u,v,feigv] = KRONREAL(a) returns the same as the previous syntax, and
% additionally returns the number of poles (finite eigenvalues of A-sE) of
% b
%
% [...] = KRONREAL(a,zerotol) returns the same as the previous syntaxes
% and uses zerotol (default: sqrt(eps)) as tolerance for 0, that is: if 
% abs(n) < zerotol, n is considered as an exact zero
%
% See also HINFCD.DSS.MINREAL. 

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

    import hinfcd.util.*;
    if nargin<2
        zerotol = sqrt(eps);
    else
        assert(isnumeric(zerotol) && isscalar(zerotol) && zerotol>=0, 'KRONREAL: zerotol must be a positive number.'); 
    end
    
    % 1. Triangular form
    % 1.1 Decomposition
    [At1,Et1,Q,Z] = qz(obj.A,obj.E,'real');
    select = (abs(ordeig(At1,Et1))<1/zerotol); feigv = sum(select); ieigv = length(Et1)-feigv;
    [At,Et,Q,Z] = ordqz(At1,Et1,Q,Z,select);
    E1 = Et(1:feigv,1:feigv);
    E2 = Et(1:feigv,(feigv+1):end);
    E3 = Et((feigv+1):end,(feigv+1):end);
    J1 = At(1:feigv,1:feigv);
    J2 = At(1:feigv,(feigv+1):end);
    J3 = At((feigv+1):end,(feigv)+1:end);
    % 1.2. Improve condition (numerical zeros -> exact zeros)
    E3(1:length(E3)+1:end) = 0; % infinite eigenvalues
    E3 = roundsmall(E3,zerotol);
    J1 = roundsmall(J1,zerotol);
    E1 = roundsmall(E1,zerotol);
    J3 = roundsmall(J3,zerotol);

    % 2. Block-diagonal form
    % 2.1 Decomposition
    [R,L] = gensylvester(E1,E3,-E2,J1,J3,-J2);
    u = blkdiag(E1,J3)\[eye(feigv) L ; zeros(ieigv,feigv) eye(ieigv)]*Q;
    v = Z*[eye(feigv) R ; zeros(ieigv,feigv) eye(ieigv)];
    E = blkdiag(eye(feigv),J3\E3);
    A = blkdiag(E1\J1,eye(ieigv));
    B = u*obj.B;
    C = obj.C*v;
    % 2.2 Improve condition (numerical zeros -> exact zeros)
    B = roundsmall(B,zerotol);
    C = roundsmall(C,zerotol); 
    
    % 3. Return object
    obj.set('A',A,'B',B,'C',C,'E',E);
    
end