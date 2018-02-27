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

function [A,B,C,D,E,Tl,Tr,K] = obsvf_ip(A,B,C,D,E,varargin)
%OBSVF_IP Transform improper system in observable form
%   Transform the system in the observability form. Works for both proper
%   and improper systems.

[A,B,C,E,Tl,Tr] = diagE_ip(A,B,C,E);
n = size(A,1);
m = size(C,1);
l = size(B,2);

%% Take care of improper states
i = 0; j = 0;
P = [A B;C D]';
Cswap = 1:m;
while(all(E(n-i,:)==0))
    if(~all(P(n-i,n+(1:m))==0))
        if(P(n-i,n+m-j)==0)
            swap = find(P(n-i,n+(1:(m-j)))~=0,1);
            if(~isempty(swap))
                P(:,[(n+swap) (n+m-j)]) = P(:,[(n+m-j) (n+swap)]);
                Cswap([(swap) (m-j)]) = Cswap([(m-j) (swap)]);
            end
        end
        for(k=1:(n-i-1))
            if(P(n-i,n+m-j)~=0)
                P(k,:) = P(k,:) - (P(k,n+m-j)*P(n-i,:)/P(n-i,n+m-j));
            end
        end
        j = j+1;
    end
    i = i+1;
end
P = P';
A = P(1:n,1:n); 
B = P(1:n,n+(1:l));
C = P(n+(1:m),1:n);
D = P(n+(1:m),n+(1:l));

%% Staircase form
if(nargin>5)
    [A,B,C,T,K] = obsvf(A,B,C,varargin{1});
else
    [A,B,C,T,K] = obsvf(A,B,C);
end
E = T*E*(T');
C = C(Cswap,:);
D = D(Cswap,:);

% TODO: Include swaps in transform
Tl = T*Tl;
Tr = Tr*(T');
end

