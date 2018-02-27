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

function umodel = USSmod(A,B,C,D)
%UTFMOD Create LFTmod from an uncertain transfer function
%   Detailed explanation goes here

    p = inputParser;
    addRequired(p,'A',@(x) assert(isnumeric(x) || isa(x,'umat') || isa(x,'ureal'),'Your matrices should be numeric or of type ''umat'' or ''ureal''.'));
    addRequired(p,'B',@(x) assert(isnumeric(x) || isa(x,'umat') || isa(x,'ureal'),'Your matrices should be numeric or of type ''umat'' or ''ureal''.'));
    addRequired(p,'C',@(x) assert(isnumeric(x) || isa(x,'umat') || isa(x,'ureal'),'Your matrices should be numeric or of type ''umat'' or ''ureal''.'));
    addRequired(p,'D',@(x) assert(isnumeric(x) || isa(x,'umat') || isa(x,'ureal'),'Your matrices should be numeric or of type ''umat'' or ''ureal''.'));
    parse(p,A,B,C,D);
       
    A = p.Results.A;
    B = p.Results.B;
    C = p.Results.C;
    D = p.Results.D;
    
    ussmodel = uss(A,B,C,D);
    umodel = UModel(ussmodel);
    
end