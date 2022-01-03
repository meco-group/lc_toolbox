function [B,V,nidxone] = niljordan(A)
% NILJORDAN Returns an upper left block Jordan form of a nilpotent matrix
%
% B = NILJORDAN(A) swaps the Jordan blocks of the nilpotent matrix A
% such that B only has a nonzero upper left block
%
% [B,V] = NILJORDAN(A) returns the same as the previous syntax, but
% additionally returns the matrix V such that B = V\A*V
%
% [B,V,nidxone] = NILJORDAN(A) returns the same as the previous syntax, but
% also returns the number of Jordan blocks with index 1

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

    assert(size(A,1)==size(A,2),'NILJORDAN: A is not square.');
    assert(all(eig(A)==0), 'NILJORDAN: A is not nilpotent.');
    n = size(A,1);
    
    if n==0
        V = zeros(0,0);
        B = V; nidxone = 0;
        return;
    end
    
    % calculate Jordan form
    [Vj,Aj] = jordan(A);
    
    % detect Jordan blocks
    d = diag(Aj,1);
    jblksiz = [];
    s = 1;
    for i=1:n-1
        if d(i)==1
            s = s+1;
        else
            jblksiz = [jblksiz s];
            s = 1;
        end
    end
    if d(end)~=1
        jblksiz = [jblksiz 1];
    else
        jblksiz = [jblksiz s];
    end
    
    % construct transformation matrix to bring index-1 blocks to the 
    % bottom right corner
    V = zeros(n);
    nidxone = sum(jblksiz==1);
    startrow = 1;
    startcol = 1;
    revstartcolidx1 = 1;
    i = 1;
    while startrow <= n
        if jblksiz(i)==1
            V(startrow,n-nidxone+revstartcolidx1) = 1;
            revstartcolidx1 = revstartcolidx1+1;
            startrow = startrow+1;
        else
            V(startrow:(startrow+jblksiz(i)-1),startcol:(startcol+jblksiz(i)-1)) = eye(jblksiz(i));
            startcol = startcol+jblksiz(i);
            startrow = startrow+jblksiz(i);
        end
        i = i+1;
    end
   
    % do the transformation
    B = V\Aj*V;
    
    % return the full transformation (might be ill-conditioned!)
    V = Vj*V;
    
    
end