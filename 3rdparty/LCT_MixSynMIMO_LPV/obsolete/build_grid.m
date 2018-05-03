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

function [grid,args] = build_grid(param,N)

if nargin == 1
    N = 5*ones(1,length(param)); % default grid
end
if isscalar(N)
    N = N*ones(1,length(param)); % same #grid points in each coordinate
end
for k = 1:length(param)
    grid{k} = linspace(param{k}.basis.domain.min,param{k}.basis.domain.max,N(k));
    args{k} = param{k}.tensor_basis.arguments;
end

end