function BAI = preim(A,S)
% PREIM Returns a basis of a preimage under a linear transformation
% Returns an orthonormal basis of the preimage of (the image of) S under 
% the linear mapping A. See S. Trenn, "Distributional differential 
% algebraic equations", PhD thesis, Institute for Mathematics, University 
% of Ilmenau, Germany, 2009.  
% 
% Helper function to make other functions more readible. 
%
% BAI = PREIM(A,S) returns BAI such that Col(BAI) spans A^(-1)(S) and BAI
% has orthonormal columns

% This file is part of hinfcd.
% Copyright (c) 2019, Laurens Jacobs, MECO Research Team @ KU Leuven. 
% 
% hinfcd is free software: you can redistribute it and/or modify it under 
% the terms of the GNU Lesser General Public License as published by the 
% Free Software Foundation, version 3.
% 
% hinfcd is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for
% more details.
% 
% You should have received a copy of the GNU Lesser General Public License
% along with hinfcd. If not, see <https://www.gnu.org/licenses/>.

    assert(size(A,1)==size(S,1), 'PREIM: A and S must have compatible dimensions.');
    BAI = null([A,S]);
    BAI = BAI(1:size(A,2),:);
    
end