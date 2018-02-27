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

function [sys,Ts,Te] = diag2(sys)
%DIAG Diagonalizes E (or A)
%   Returns an equivalent system with E (or A) equal to a (semi)-unity
%   matrix.
%
%   INPUTS:
%   sys = system to be transformed
%
%   OUTPUTS:
%   sys: the transformed system
%   Ts,Te: Enew = Te*Eold*Ts, Anew = Te*Aold*Ts, Bnew = Te*Bold, Cnew =
%   Cold*Ts

Ts = eye(size(sys.a));
Te = eye(size(sys.a));

if ~isempty(sys.e)
    [U,S,V] = svd(sys.e);
    scal = diag(1./diag(S));
    scal(~isfinite(scal)) = 1;

    sys.e = scal*transpose(U)*sys.e*V;
    sys.a = scal*transpose(U)*sys.a*V;
    sys.b = scal*transpose(U)*sys.b;
    sys.c = sys.c*V;
    
    Ts = V;
    Te = transpose(U);
else
    sys.e = eye(size(sys.a));
end

end

