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

function [sys] = sstransform(sys,L,R)
%SSTRANSFORM
%   Transform the system using left transform matrix L and right
%   transform R. The system will be of the form:
%   L*E*R*dot(x) = L*A*R*x + L*B*u;
%   y = C*R*x + D*u;

if((all(size(L)==size(R)))&&(all(size(L)==size(sys.a))))
    sys.e = L*sys.e*R;
    sys.a = L*sys.a*R;
    sys.b = L*sys.b;
    sys.c = sys.c*R;
else
    error('L and/or R are of incompatible size.');
end

end

