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

clear all
close all
clc

%% 1. Goal
% The goal of this tutorial is to design a controller based on an
% H-infinity formalism. Depending on the type of constraints we impose, a
% different solver is most suited and will therefore be picked by the
% toolbox.

%% 2. System declarations
% Plant
Gmod = ZPKmod([-16+720j, -16-720j],[-1100, 190, -160, -18+770j, -18-770j],1e8);

% Weights
MS = Weight.DC(10);          % Maximum on sensitivity: 8dB
WS = Weight.LF(20,1,-40);             % Weight on sensitivity to assure roll-off
WSu = Weight.LF(20,1);
WU = Weight.HF(2e3,1,-10);   % Weight on input sensitivity to assure roll-off on the controller
WD = Weight.DC(10);                % Problem with constraint does not solve :(

%% 3. Controller design
G = IOSystem(1,1);
G.add(Gmod);
K = IOSystem(1,1);

% Declare the external reference so we can use it
r = Signal();
k = Signal();
d = Signal();
    
% assign some convenient variables
c = k + d;
u = G.in;
y = G.out;
e = r - y;
connections = [G.in == c; K.in == e; K.out == k];
P = IOSystem(G,K,connections);
    
% MIMO design: take disturbance rejection into account
S = Channel(e/r, 'Sensitivity');
U = Channel(u/r, 'Input Sensitivity');
D = Channel(y/d, 'Disturbance Rejection');
    
objective = [WS*S];
constraints = [WU*U <= 1, MS*S <= 1];
options = struct('controller_name','mixed\_controller');
[P,C1,info1] = P.solve(objective, constraints, K, options);

objective = [WS*S + WD*D];
constraints = [WU*U <= 1, MS*S <= 1];
options = struct('FullOrderSolver','mixedHinfsynMIMO','controller_name','mixed\_controller\_mimo-dist_rej');
[P,C2,info2] = P.solve(objective, constraints, K, options);

bodemag(info1,info2,S,U,D)

%% 4. Discussion
% Designing a controller like this, we can decide also on the disturbance
% rejection. This however involves a new input 'd' which represents a
% disturbance in our control loop.
% 
% Because of the second exogeneous input, the controller requires a MIMO
% solver, which is also implemented in mixedHinfsynMIMO.
%
% The resulting problem is badly conditioned as can be seen from the
% comparison of the weight on U and the resulting transfer function.