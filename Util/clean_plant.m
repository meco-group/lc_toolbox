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

function [sys] = clean_plant(sys,varargin)
%CLEAN_PLANT Summary of this function goes here
%   Detailed explanation goes here

if(nargin>1)
    tol = varargin{1};
else
    tol = eps;
end

members = {'a','b','c','d','e'};
for i = 1:length(members)
    field = sys.(members{i});
    field(abs(field) < tol) = 0;
    sys.(members{i}) = field;
end

end

