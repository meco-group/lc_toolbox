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

function poly = butterpoly(n,u)
%BUTTERPOLY Make butterworth polynomial
%   Constructs the coefficients of a butterworth polynomial of degree n in
%   continuous time. 
%   Args: n = degree
%         u = abs(roots)

if mod(n,2)==1 %uneven
    poly = [1 u];
    for k = 1:((n-1)/2)
        poly = conv(poly,[1 -2*u*cos((2*k+n-1)*pi/(2*n)) u^2]);
    end
else
    poly = 1;
    for k = 1:(n/2)
        poly = conv(poly,[1 -2*u*cos((2*k+n-1)*pi/(2*n)) u^2]);
    end
end

end