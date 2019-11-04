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

function umodel = UTFmod(num, den)
%UTFMOD Create LFTmod from an uncertain transfer function
%   Detailed explanation goes here

    p = inputParser;
    addRequired(p,'num',@(x) assert(isnumeric(x) || isa(x,'umat') || isa(x,'ureal'),'Your numerator should be numeric or of type ''umat'' or ''ureal''.'));
    addRequired(p,'den',@(x) assert(isnumeric(x) || isa(x,'umat') || isa(x,'ureal'),'Your denominator should be numeric or of type ''umat'' or ''ureal''.'));
    parse(p,num,den);
       
    num = p.Results.num;
    den = p.Results.den;
    
    ussmodel = uss(tf(num,den));
    umodel = UModel(ussmodel);
    
end