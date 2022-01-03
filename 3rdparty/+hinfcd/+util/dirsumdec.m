function Z = dirsumdec(varargin)
% DIRSUMDEC Returns the summand required to decompose a vector space through the direct sum
%
% Z = DIRSUMDEC(A,B,C,...,sum) returns Z, with orthonormal columns, such 
% that the direct sum of Col(A), Col(B), ... and the subspace Col(Z)
% forms Col(sum). 
%
% Helper function to make other functions more readible. 

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

    import hinfcd.util.*;

    assert(nargin>=2,'DIRSUMDEC: Not enough input arguments.');
    assert(all(cellfun(@(x) isnumeric(x) && size(x,1)==size(varargin{1},1), varargin)),'DIRSUMDEC: Inputs must be matrices with the same numbers of rows.');
    sum = varargin{end};
    
    % check for linear dependency 
    m = horzcat(varargin{1:end-1});
    assert(size(m,2)<=size(m,1) && rank(m')<=size(m,2),'DIRSUMDEC: The summands provided are not linearly independent. A direct sum decomposition does not exist.');
        
    % calculate the summand that makes Z = R^n (n = size(sum,1))
    Z = null(m');
    
    % remove the dimensions of R^n such that the columns of Z span exactly
    % the same space as the columns of sum
    Z = isectspace(Z,sum);

end