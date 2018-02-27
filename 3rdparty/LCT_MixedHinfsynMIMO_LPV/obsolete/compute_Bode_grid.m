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

function compute_Bode_grid(sys,param,N)

[grid,args] = build_grid(param,N);      % parameter grid
index_list  = build_indices_grid(grid); % index list grid
M = sys.M.grid_eval(grid,args);         % gridded system matrix
figure; % draw bode plots
for k = 1:size(index_list,1)
    bode(LFTmod(squeeze(M(index_list{k,:},:,:)),sys.Nu,sys.Nl,eye(sys.nx))); hold on
end

end