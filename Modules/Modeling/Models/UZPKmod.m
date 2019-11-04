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

function umodel = UZPKmod(z,p,k)
%UZPKMOD Create LFTmod from an uncertain zpk model
%   Detailed explanation goes here

    pr = inputParser;
    addRequired(pr,'z',@(x) assert(isnumeric(x) || isa(x,'umat') || isa(x,'ureal'),'Your zeros should be numeric or of type ''umat'' or ''ureal''.'));
    addRequired(pr,'p',@(x) assert(isnumeric(x) || isa(x,'umat') || isa(x,'ureal'),'Your poles should be numeric or of type ''umat'' or ''ureal''.'));
    addRequired(pr,'k',@(x) assert((isnumeric(x) || isa(x,'umat') || isa(x,'ureal')) && length(x) == 1,'Your gain should be a scalar, either numeric or of type ''umat'' or ''ureal''.'));
    parse(pr,z,p,k);
       
    z = pr.Results.z;
    p = pr.Results.p;
    k = pr.Results.k;
    
    s = tf('s');    
    num = 1; den = 1;
    
    for i = 1:length(z)
        num = num*(s-z(i));
    end
    for i = 1:length(p)
        den = den*(s-p(i));
    end
    
    ussmodel = uss(k*num/den);
    umodel = UModel(ussmodel);
    
end