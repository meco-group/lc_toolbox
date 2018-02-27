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

function [sysout] = ssminreal(sysin,varargin)
%SSMINREAL Transform a state space model into a minimal description
%   SSMINREAL transforms a general state space model (proper and improper)
%   into a minimal description. If an improper system can be represented as with a
%   proper state space model, SSMINREAL will transform it into a proper
%   form.
%   input:
%   sys = (descriptor) state space model to be converted
%   output:
%   sys = minimal representation of the original model
%   Tl, Tr = transformation matrices to retrieve the staircase
%   representation of the model, from which a minimal representation can be
%   deduced:
%   At = Tl*A*Tr; Bt = Tl*B;
%   Ct = C*Tr; Dt = D;
%   Et = Tl*E*Tr;
%   If the input model is of the descriptor form, but the underlying system
%   is proper, the transformation matrices are not correct!
%   TODO: Add correct transformation from dss2ss;

if nargin > 1
    tol = varargin{1};
else
    tol = 1e-4;
end

flag = false;
if isa(sysin,'AbstractDSSmod')
    sys = std(sysin);
    flag = true;
else
    sys = sysin;
end

[sys,is_ss,rtot] = dss2ss(sys); % Transform a possibly improper description to a possibly proper description

A = sys.A; B = sys.B;
C = sys.C; D = sys.D;
E = sys.E;
n = size(A,1);

Tl = eye(size(A)); Tr = eye(size(A)); 

%% Reduce the system to an observable form
[A,B,C,D,E,Tlo,Tro,~] = obsvf_ip(A,B,C,D,E);
P = [A B;C D];
r = rmin(P',n,tol); 
rtot = rtot+r;

Tl = Tlo*Tl; Tr = Tr*Tro;

P(1:r,:) = []; P(:,1:r) = []; E(1:r,:) = []; E(:,1:r) = [];
n = n-r;
A = P(1:n,1:n); B = P(1:n,n+1:end);
C = P(n+1:end,1:n); D = P(n+1:end,n+1:end);

%% Reduce the system to a controlable form
[A,B,C,D,E,Tlc,Trc,~] = ctrbf_ip(A,B,C,D,E); 
P = [A B;C D];
r = rmin(P,n,tol);
rtot = rtot+r;

Tl((end-n+1):end,:) = Tlc*Tl((end-n+1):end,:); Tr(:,(end-n+1):end) = Tr(:,(end-n+1):end)*Trc;

P(1:r,:) = []; P(:,1:r) = []; E(1:r,:) = []; E(:,1:r) = [];
n = n-r;
A = thr(P(1:n,1:n)); B = thr(P(1:n,n+1:end));
C = thr(P(n+1:end,1:n)); D = thr(P(n+1:end,n+1:end));
E = thr(E);

%% Extract the system
if(is_ss)
    disp('ssminreal: The output system is proper.');
else
    disp('ssminreal: The output system is improper.');
end
if rtot > 0
    fprintf('ssminreal: %d state(s) removed.\n',rtot);
    if flag
        sysout = sysin;
        sysout = sysout.setdssdata(A,B,C,D,E);
    else
        sysout = dss(A,B,C,D,E,sysin);
    end
else
    sysout = sysin;
end
end

function r = rmin(P,n,tol)
r = 0;
for k=1:n
    P0 = P(1:k,(k+1):end);
    P0 = P0(:);
    if(all(abs(P0)<=tol))
        r = k;
    end
end
end

function m = thr(m,tol)
    if nargin == 1
        tol = eps;
    end
    m(abs(m) < tol) = 0;
end
