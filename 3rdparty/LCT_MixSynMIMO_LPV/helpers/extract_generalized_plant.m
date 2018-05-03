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

function gen_sys = extract_generalized_plant (sys,nu,ny)
% extract the generalized plant

            if nargin == 1; gen_sys.nu = 0; gen_sys.ny = 0; end
            if nargin == 2; gen_sys.ny = 0; end
            gen_sys.A = sys.A;
            gen_sys.nx = size(gen_sys.A,1);
            gen_sys.nu = nu;
            gen_sys.nw = size(sys.B,2)-nu;
            gen_sys.ny = ny;
            gen_sys.nz = size(sys.C,1)-ny;
            if gen_sys.nu ~= 0 && gen_sys.ny == 0     % state feedback case
                gen_sys.Bw = temporary_fix_subsref(sys.B,1:gen_sys.nx,1:gen_sys.nw); gen_sys.Bu = temporary_fix_subsref(sys.B,1:gen_sys.nx,gen_sys.nw+1:gen_sys.nw+gen_sys.nu);
                gen_sys.Dw = temporary_fix_subsref(sys.D,1:gen_sys.nz,1:gen_sys.nw); gen_sys.Du = temporary_fix_subsref(sys.D,1:gen_sys.nz,gen_sys.nw+1:gen_sys.nw+gen_sys.nu);
            elseif gen_sys.nu ~= 0 && gen_sys.ny ~= 0 % output feedback case
                gen_sys.Bw  = temporary_fix_subsref(sys.B,1:gen_sys.nx,1:gen_sys.nw);
                gen_sys.Bu  = temporary_fix_subsref(sys.B,1:gen_sys.nx,gen_sys.nw+1:gen_sys.nw+gen_sys.nu);
                gen_sys.Cz  = temporary_fix_subsref(sys.C,1:gen_sys.nz,1:gen_sys.nx);
                gen_sys.Dzw = temporary_fix_subsref(sys.D,1:gen_sys.nz,1:gen_sys.nw);
                gen_sys.Dzu = temporary_fix_subsref(sys.D,1:gen_sys.nz,gen_sys.nw+1:gen_sys.nw+gen_sys.nu);
                gen_sys.Cy  = temporary_fix_subsref(sys.C,gen_sys.nz+1:gen_sys.nz+gen_sys.ny,1:gen_sys.nx);
                gen_sys.Dyw = temporary_fix_subsref(sys.D,gen_sys.nz+1:gen_sys.nz+gen_sys.ny,1:gen_sys.nw);
                gen_sys.Dyu = temporary_fix_subsref(sys.D,gen_sys.nz+1:gen_sys.nz+gen_sys.ny,gen_sys.nw+1:gen_sys.nw+gen_sys.nu);
            end
            
            function part = temporary_fix_subsref(M,rows,cols)
                % construct first row
                part = M(rows(1),cols(1)); % upper-left element
                for i = 2:length(cols) % concatenate columns to construct 1st row
                    part = [part, M(rows(1),cols(i))];
                end
                
                % construct row j and concatenate to previously constructed rows
                for j = 2:length(rows)
                    rowj = M(rows(j),cols(1)); % upper-left element
                    for i = 2:length(cols) % concatenate columns to construct 1st row
                        rowj = [rowj, M(rows(j),cols(i))];
                    end
                    part = [part; rowj];
                end
            end
end