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

function [sys,T] = rrefE(sys)
%RREFE Put the system dynamics in a canonical form
%   Returns an equivalent system with E in a reduced row echelon form
%
%   INPUTS:
%   sys = system to be transformed
%
%   OUTPUTS:
%   sys: the transformed system
%   T: Enew = T*Eold, Anew = T*Aold, Bnew = T*Bold


[A,B,C,D,E] = dssdata(sys);
T = eye(size(A));

in = [E A B T];
out = rref(in);

e = size(E,2);
b = size(B,2);

sys.E = out(:,1:e);
sys.A = out(:,e+1:2*e);
sys.B = out(:,2*e+1:2*e+b);
T = out(:,end-e+1:end);

end

