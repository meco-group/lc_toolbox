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
        if isnumeric(u)
            poly = conv(poly,[1 -2*u*cos((2*k+n-1)*pi/(2*n)) u^2]);
        else
            poly = convpar(poly,[1 -2*u*cos((2*k+n-1)*pi/(2*n)) u^2]);
        end
    end
else
    poly = 1;
    for k = 1:(n/2)
        if isnumeric(u)
            poly = conv(poly,[1 -2*u*cos((2*k+n-1)*pi/(2*n)) u^2]);
        else
            poly = convpar(poly,[1 -2*u*cos((2*k+n-1)*pi/(2*n)) u^2]);   
        end
    end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [y] = convpar(u,v)
check_optispline;
import splines.*;
opti = OptiSplineYalmip();
TB = TensorBasis({},{});

nu = size(u,2);
nv = size(v,2);
n = nu+nv-1;
%y = {};
% y = zeros(1,n);
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
U = [u zeros(1,n-nu)];
V = [v zeros(1,n-nv)];
ys = {};
for i = 1:n
        x = 0;
    for j = 1:i
        x = x + U(j)*V(i-j+1);
    end
    ys = {ys{:}, x};
end
    y = [ys{:}];
end



