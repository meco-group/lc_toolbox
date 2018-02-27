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

function [ dPdt ] = build_derivative(P,varargin)
% This function computes the derivative of any B-spline function and
% returns zeros if the variable is constant and returns rate dependent
% function in case of varying variable.
n = size(P,1);
dPdt = zeros(n);
if (nargin == 2)
    param = varargin{1};
for i = 1:length(param)
    if ~(all(param{i}.rate == 0)) && ~any(isinf(param{i}.rate)) % if rate is not zero or unbounded, take derivative
        p = SchedulingParameter(strcat(param{i}.tensor_basis.arguments,'dot'),param{i}.rate,0);
        dPdt = dPdt + P.derivative(1,param{i}.tensor_basis.arguments)*p;
    end
end
end
end

