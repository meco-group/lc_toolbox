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

clear all
close all
clc

A = [-100 0 0;50000000 -100 100000;100 0 -0.002];
fprintf('Original condition number of A: %e (rcond = %e)\n',cond(A),rcond(A));
At = A;

for k=1:40
As = (At+At')/2;

[P,D] = eig(As); Dq = diag(diag(D).^(1/4));

T = P*Dq*(P'); Tinv = P*pinv(Dq)*(P');
At = Tinv*At*T;

fprintf('Current condition number of A: %e (rcond = %e)\n',cond(At),rcond(At));
end