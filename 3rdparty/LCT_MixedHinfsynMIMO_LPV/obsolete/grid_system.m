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

function grid_sys = grid_system(sys,param,N)

[grid,args] = build_grid(param,N);      % parameter grid
index_list  = build_indices_grid(grid); % index list grid

% grid system
if isa(sys.A,'splines.Function')
    grid_sys.A = sys.A.grid_eval(grid,args);
else
    for k = 1:size(index_list,1)
        grid_sys.A(index_list{k,:},:,:) = sys.A;
    end
end
if isa(sys.B,'splines.Function')
    grid_sys.B = sys.B.grid_eval(grid,args);
else
    for k = 1:size(index_list,1)
        grid_sys.B(index_list{k,:},:,:) = sys.B;
    end
end
if isa(sys.C,'splines.Function')
    grid_sys.C = sys.C.grid_eval(grid,args);
else
    for k = 1:size(index_list,1)
        grid_sys.C(index_list{k,:},:,:) = sys.C;
    end
end
if isa(sys.D,'splines.Function')
    grid_sys.D = sys.D.grid_eval(grid,args);
else
    for k = 1:size(index_list,1)
        grid_sys.D(index_list{k,:},:,:) = sys.D;
    end
end

end
