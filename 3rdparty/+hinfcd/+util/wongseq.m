function W = wongseq(E,A,B,C,S)
% WONGSEQ Returns the limit of a Wong sequence
% Returns the Wong sequence that is required to calculate the Kalman
% decomposition of a descriptor system using its geometric properties.
% 
% The Wong sequence of subspace Col(S) we consider is recursively defined as:
%   W(0) = S
%   W(i+1) = intersect(preim(A,E*W(i)) + im(B)), ker(C))
% where Col(W(n)), n = length(E), is the subspace of interest.
%
% The approach is adopted from A. Banaszuk, M. Kociecki and F.L. Lewis,
% "Kalman Decomposition for Implicit Linear Systems", IEEE Transactions on
% Automatic Control, vol. 37, no. 10, October 1992. 
%
% W = WONGSEQ(E,A,B,C,S) returns the converged Wong sequence, i.e. W 
% satisfies W = WONGSEQ(E,A,B,C,W)

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

    assert(size(E,1)==size(E,2) && size(A,1)==size(A,2) && size(A,1)==size(E,1), 'WONGSEQ: E and A must be square matrices of the same size.');
    assert(size(A,1)==size(B,1), 'WONGSEQ: A and B should have the same number of rows.'); 
    assert(size(A,2)==size(C,2), 'WONGSEQ: A and C should have the same number of columns.'); 
    assert(size(A,1)==size(S,1), 'WONGSEQ: S should have the same number of rows as A and B.'); 
    
    % calculate constants
    imB = orth(B); 
    kerC = null(C);
    
    % initial element
    W = S; 

    % iterate
    for i=1:size(A,1)+1
        Wn = isectspace(preim(A,orth([E*W,imB])),kerC);
        if size(Wn,2)==size(W,2)
            break; % converged
        elseif i==size(A,1)+1
            assert(size(Wn,2)==size(W,2), 'WONGSEQ: Wong sequence did not converge.');
        end
        W = Wn;
    end
    
end