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

function [sys,Hl,Hr] = ssStaircase(sys)
%SSSTAIRCASE
%   Implementation of the staircase algorithm of M.M. Rosenbrock,
%   State Space and Multivariable Theory 

%% Step 1: Reshape E in proper form:
%   1   0   0   ..  0
%   0   1   0   ..  0   
%   0   0   1   ..  0
%   |   |   |   ..  0
%   0   0   0   ..  0
n = size(sys.a,1);
m = size(sys.c,1);
l = size(sys.b,2);



% if(~isempty(sys.e))
%     R = rref([sys.e eye(n)]);
%     sys.e = R(1:n,1:n);
%     Hl = R(:,(n+1):(2*n));
%     sys.a = Hl*sys.a;
%     sys.b = Hl*sys.b;
%     
%     R = rref([sys.e' eye(n)]);
%     sys.e = R(1:n,1:n);
%     Hr = R(:,(n+1):(2*n))';
%     sys.a = sys.a*Hr;
%     sys.c = sys.c*Hr;
%     zpk(sys)
% else
%     Hl = eye(n);
%     Hr = eye(n);
% end

%% Step 2: Put A in staircase form
% i = 0; j = 0;
% eps = 1e-12;
% m = 1;
% 
% while(((i-j)<l)&&(j~=n)) % Stop condition: all rows are reduced or finished
%     sys = cleanSys(sys);
%     P = sys2mat(sys);
%     if(all(abs(P(1:(n-j),n+l-i))<eps)) % all elements are zero
%         i = i+1; % step 6
%     else
%         if(P(n-j,n+l-i)==0) % Check if the last element in the column is nonzero 
%             ind = find(abs(P(1:n-j,n+l-i))>eps,1); % Find a nonzero element in the column to swap the last element with
%             [sys,Hi,Hiinv] = swap(sys,ind,n-j); 
%             P = sys2mat(sys);
%             Hr = Hr*Hi; Hl = Hiinv*Hl; % Update transformation matrices
%         end
%         for(k=1:(n-j-1)) % Reduce all other elements to zero
%             if(abs(P(k,n+l-i))>eps)
%                 [sys,Hi,Hiinv] = add(sys,k,n-j,-P(k,n+l-i)/P(n-j,n+l-i));
%                 P = sys2mat(sys);
%                 Hr = Hr*Hi; Hl = Hiinv*Hl;
%             end
%         end
%         
%         i = i+1; % step 5
%         j = j+1;
%     end
% end

[Ao,Bo,Co,Eo,Tlo,Tro,Ko] = obsvf_ip(sys.A,sys.B,sys.C,sys.E);
[Ac,Bc,Cc,Ec,Tlc,Trc,Kc] = ctrbf_ip(sys.A,sys.B,sys.C,sys.E);

eps = 1e-6;
% [Ao,Bo,Co,To] = obsvf(sys.a,sys.b,sys.c);
% [Ac,Bc,Cc,Tc] = ctrbf(sys.a,sys.b,sys.c);
Po = [Ao Bo;Co sys.d]; Po = Po';
Pc = [Ac Bc;Cc sys.d]
Po(abs(Po)<eps) = 0
% Pc = P 
Pc(abs(Pc)<eps) = 0;

%% Go to minimal description
ro = minimise(Po,n)
rc = minimise(Pc,n)
if(ro>rc)
    r = ro
    sys.e = Eo;
    Hl = Tlo; Hr = Tro;
    P = Po';
else
    r = rc
    sys.e = Ec;
    Hl = Tlc; Hr = Trc;
    P = Pc;
end

P(1:r,:) = []; P(:,1:r) = [];
n = n-r;
sys = dss(P(1:n,1:n),P(1:n,n+1:end),P(n+1:end,1:n),P(n+1:end,n+1:end),sys.e(r+1:end,r+1:end));

end

function [P] = sys2mat(sys)
% Compose system matrix (A B;C D)

P = zeros(size(sys.a,1)+size(sys.c,1),size(sys.a,2)+size(sys.b,2));
P(1:size(sys.a,1),1:size(sys.a,2)) = sys.a;
P(size(sys.a,1)+[1:size(sys.c,1)],1:size(sys.a,2)) = sys.c;
P(1:size(sys.a,1),size(sys.a,2)+[1:size(sys.b,2)]) = sys.b;
P(size(sys.a,1)+[1:size(sys.c,1)],size(sys.a,2)+[1:size(sys.b,2)]) = sys.d;
end

function [sys] = cleanSys(sys)
eps = 1e-6;
sys.a(abs(sys.a)<eps) = 0;
sys.b(abs(sys.b)<eps) = 0;
sys.c(abs(sys.c)<eps) = 0;
sys.d(abs(sys.d)<eps) = 0;
sys.e(abs(sys.e)<eps) = 0;
end

function [rf] = minimise(P,n)
rf = 0;
for r=1:n
    P0 = P(1:r,(r+1):end);
    P0 = P0(:);
    if(all(abs(P0)<=1e-9))
        rf = r;
    end
end
end

function [sys,H,Hinv] = swap(sys,p,q)
n = size(sys.a,1);
m = size(sys.c,1);
l = size(sys.b,2);

if((p>n)||(q>n))
    error('Swap indices are out of bounds!');
end
if(p>q)
    t = p;
    p = q;
    q = t;
elseif(p==q)
    error('Swap indices must not be equal!');
end

H = eye(n);
Hinv = eye(n);
Hs = eye(n+l);
Hsinv = eye(n+m);

Hsub = eye(q-p+1);
Hsub(1,1) = 0; Hsub(1,q-p+1) = 1;
Hsub(q-p+1,1) = 1; Hsub(q-p+1,q-p+1) = 0; 
Hsubinv = Hsub; % Check!
% Hsubinv*Hsub

H(p:q,p:q) = Hsub;
Hinv(p:q,p:q) = Hsubinv;

sys.a = Hinv*sys.a*H;
sys.b = Hinv*sys.b;
sys.c = sys.c*H;
if(~isempty(sys.e))
    sys.e = Hinv*sys.e*H;
end

end

function [sys,H,Hinv] = add(sys,p,q,gamma)

n = size(sys.a,1);
m = size(sys.c,1);
l = size(sys.b,2);

H = eye(n,n); H(p,q) = -gamma;
Hinv = eye(n,n); Hinv(p,q) = gamma;

sys.a = Hinv*sys.a*H;
sys.b = Hinv*sys.b;
sys.c = sys.c*H;
if(~isempty(sys.e))
    sys.e = Hinv*sys.e*H;
end

end