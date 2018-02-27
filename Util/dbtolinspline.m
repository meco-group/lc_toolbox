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

function spline_out = dbtolinspline(spline_in,noiknots,degree)
%DBTOLINSPLINE Approximates 10^(spline_in/20) as a new B-spline in the same
%arguments as spline_in
%   Detailed explanation goes here

    import splines.*;

    if nargin == 1
        noiknots = 3;
        degree = 3;
    elseif nargin == 3
    else error('dbtolinspline needs one or three input arguments.');
    end
    
    % create a new tensor basis 
        
    domains = cellfun(@(x) x.domain.data, spline_in.tensor_basis.bases, 'un', 0);
    bases = cellfun(@(x) BSplineBasis(x, degree, noiknots),domains,'un',0);
    new_basis = TensorBasis(bases,cellstr(spline_in.tensor_basis.arguments)');
    
    % fit the new spline based on a linspace of 100 points in each
    % dimension
    
    grid = cellfun(@(x) linspace(x(1),x(2),100),domains,'un',0);
    args = cellstr(spline_in.tensor_basis.arguments)';
    val = spline_in.grid_eval(grid,args);
    val = 10.^(val/20);

    opti = OptiSpline();
    F = opti.Function(new_basis);
    if length(bases) == 1
        tmp = F.grid_eval(grid,args);
        diff = tmp(:,:,1) - val;
    else
        diff = squeeze(F.grid_eval(grid,args)) - val;
    end
    vars = diff.data;
    obj = norm(vars,2); % least squares
    
    opti.minimize(obj);% unconstrained: boundaries might not be the same anymore!
    opti.solver('ipopt');
    sol = opti.solve();
    F = sol.value(F);
    
    spline_out = F;
    
end

