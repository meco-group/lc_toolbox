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

function plot_bound(f,param,N)  

[grid,args] = build_grid(param,N);
if isa(f,'splines.Function')
    f = f.grid_eval(grid,args);      
end
if length(grid) == 1
    xlabel(args{1}); ylabel('bound');
    plot(grid{1},f); 
elseif length(grid) == 2
    xlabel(args{1}); ylabel(args{2}); zlabel('bound');
    surf(grid{1},grid{2},f');
else
    error('plotting works only for functions depending on 1 or 2 variables')
end

end
