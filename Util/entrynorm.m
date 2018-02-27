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

function output = entrynorm(S,varargin)
%ENTRYNORM Computes an entrywise norm
%   Computes the system norm entrywise. All usual norm options are
%   available (basicly choosing the norm: 2 or inf)

[m,n] = size(S);
output = zeros(m,n);

warning off;
for i = 1:m
    for j = 1:n
        output(i,j) = norm(S(i,j),varargin{:});
    end
end
warning on;

end

