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
clear global
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
MS = Weight.DC(10);         % Maximum on sensitivity: 8dB
WS = Weight.LF(80,1,-100);  % Weight on sensitivity to assure roll-off
MU = Weight.DC(20);         % Maximum input sensitivity: 20dB
WU = Weight.HF(1e3,1,-60);  % Weight on input sensitivity to assure roll-off on the controller
beta = 5;

%% 3. Controller design
G = IOSystem(1,1);
G.add(Gmod);
K = IOSystem(1,1);  

r = Signal();
u = G.in;
y = G.out;
e = r - y;
connections = [K.in == e; K.out == u];
P = IOSystem(G,K,connections);

% 'Classical' vs 'mixed' H-inf design
S = Channel(e/r, 'Sensitivity');
U = Channel(u/r, 'Input Sensitivity');
    
objective = [WS*S;WU*U];
constraints = [];
options = struct('controller_name','stacked\_controller');
[P,C1,info1] = P.solve(objective, constraints, K, options);

objective = [WS*S];
constraints = [WU*U <= 1];
options = struct('controller_name','mixed\_controller');
[P,C2,info2] = P.solve(objective, constraints, K, options);

bodemag(info1,info2,S,U)

%% 4. Discussion
% Switching to a mixed design eases the problem formulation. Objectives and
% constraints can be rigorously separated contrary to the 'classical' H-inf 
% design where the only specification is an objective.
%
% Choosing the right bandwidth for the sensitivity, we can come up with a
% nearly identical controller as we did for the first case. However, it is
% now easier to just play around with the desired bandwidth on the
% sensitivity function and come up with a desirable response.
%
% Maybe play around a bit with beta and try to explain what happens to the
% controller.