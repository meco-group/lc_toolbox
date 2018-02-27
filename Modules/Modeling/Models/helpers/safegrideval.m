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

function [e] = safegrideval(func,grid,args)
% Ensures safe evaluation of a model matrix on a grid. 
% Written for the Optispline toolbox. 
%
% Parameters: 
%  func : object to be evaluated, either \c double or \c Function
%  grid : grid : cell with grid{i} a vector (\c double) with the gridpoints for
%  scheduling parameter i
%  args : cell with args{i} the name (\c char) of scheduling parameter i @type cell

    grid = cellfun(@(x) x(:),grid,'un',0);

    if isnumeric(func)
        r = repmat(func,[1,1,cellfun(@length,grid)]);
        e = permute(r,[2+(1:length(grid)),1,2]);
    else
        e = func.grid_eval(grid,args);
    end
end