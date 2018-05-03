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

function grid_norm = compute_norm_grid(sys,ny,nu,param,N,opts,norm)

grid_sys   = grid_system(sys,param,N); % gridded system
grid       = build_grid(param,N);      % parameter grid
index_list = build_indices_grid(grid); % index list grid
for k = 1:size(index_list,1)
    sys_loc.A = squeeze(grid_sys.A(index_list{k,:},:,:));
    sys_loc.B = squeeze(grid_sys.B(index_list{k,:},:,:));
    sys_loc.C = squeeze(grid_sys.C(index_list{k,:},:,:));
    sys_loc.D = squeeze(grid_sys.D(index_list{k,:},:,:));
    G = extract_generalized_plant(sys_loc,nu,ny);
    if norm == inf
        Info = compute_FO_Hinf_LTI(G.A,G.Bw,G.Bu,G.Cz,G.Cy,G.Dzw,G.Dzu,G.Dyw,G.Dyu,opts); 
        grid_norm(index_list{k,:}) = Info.gam;
    elseif norm == 2
        error('gridded H2 norm not yet supported')
    else
        error('select Hinf or H2 norm');
    end
end

end