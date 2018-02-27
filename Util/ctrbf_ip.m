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

function [A,B,C,D,E,Tl,Tr,K] = ctrbf_ip(A,B,C,D,E,varargin)
%CTRBF_IP Transform improper system in controlable form
%   Transform the system in the controlability form. Works for both proper
%   and improper systems.

[A,B,C,E,Tl,Tr] = diagE_ip(A,B,C,E);
n = size(A,1);
m = size(C,1);
l = size(B,2);

sys = dss(A,B,C,D,E);

%% Take care of improper states
% i = 0; j = 0;
% P = [A B;C D];
% Bswap = (1:l);
% while(all(E(n-i,:)==0))
%     if(P(n-i,n+l-j)==0)
%         swap = find(P(n-i,n+(1:(l-j)))~=0,1);
%         if(~isempty(swap))
%             P(:,[(n+swap) (n+l-j)]) = P(:,[(n+l-j) (n+swap)]);
%             Bswap([(swap) (l-j)]) = Bswap([(l-j) (swap)]);
%         else
%             warning('no input to the algebraic equation. We have to redefine the states');
% %             replace = find(abs(P(n-i,1:n))>eps,1);
% %             replace = 4;
% %             T = eye(n); T(replace,:) = P(n-i,1:n);
% %             E = T\E*T; E(replace,replace) = 0; E(n,:)=[]; E(:,replace+1)=[];
% %             A = T\A*T; A(n,:)=[]; A(:,replace-1)=[];
% %             B = T\B; B(n,:) = [];
% %             C = C*T; C(:,replace-1) = [];
% %             sysred = dss(A,B,C,D,E);
% %             figure, bode(sys,sysred);
% %             break;
%             P(:,[(n+swap) (n+l-j)]) = P(:,[(n+l-j) (n+swap)]);
%             Bswap([(swap) (l-j)]) = Bswap([(l-j) (swap)]);
%         end
%     end
%     for(k=1:(n-i-1))
%         if(P(n-i,n+l-j)~=0)
%             P(k,:) = P(k,:) - (P(k,n+l-j)*P(n-i,:)/P(n-i,n+l-j));
%         end
%     end
%     i = i+1;
% end
% A = P(1:n,1:n); 
% B = P(1:n,n+(1:l));
% C = P(n+(1:m),1:n);
% D = P(n+(1:m),n+(1:l));
Bswap = 1:l;

%% Staircase form
if(nargin>5)
    [A,B,C,T,K] = ctrbf(A,B,C,varargin{1});
else
    [A,B,C,T,K] = ctrbf(A,B,C);
end
E = T*E*(T');
B = B(:,Bswap);
D = D(:,Bswap);

Tl = T*Tl;
Tr = Tr*(T');
end

