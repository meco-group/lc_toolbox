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

function [A,B,C,E,Tl,Tr] = diagE_ip(A,B,C,E)
%DIAGE_IP
%   Put plant in diagonal form. This results in a matrix of the form: 
%   E = [I 0;0 0]. Straightforward use can be diagonalizing E, but it can also
%   be used to diagonalize A, by switching inputs A and E.
%   inputs: system matrices A,B,C,E
%   outputs: transformed At,Bt,Ct,Et where E is in the diagonal form and the
%   transformation matrices Tl and Tr so that Tl*A*Tr = At, Tl*B = Bt, C*Tr
%   = Ct, Tl*E*Tr = Et

% n = size(A,1);
% 
% if(~isempty(E))
%     R = robust_rref([E eye(n)]);
%     E = R(1:n,1:n);
%     Tl = R(:,(n+1):(2*n));
%     A = Tl*A;
%     B = Tl*B;
%     
%     R = robust_rref([E' eye(n)]);
%     E = R(1:n,1:n);
%     Tr = R(:,(n+1):(2*n))';
%     A = A*Tr;
%     C = C*Tr;
% else    
%     E = eye(n);
%     Tl = eye(n);
%     Tr = eye(n);
% end

sys = dss(A,B,C,zeros(size(C,1),size(B,2)),E);
[sys,Ts,Te] = diag2(sys);

A = sys.a;
B = sys.b;
C = sys.c;
E = sys.e;
Tl = Te;
Tr = Ts;

end

