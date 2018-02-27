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
figure('Name','Plant bode diagram');
bode(Gmod);

% Weights
MS = Weight.DC(-10);        % Maximum on sensitivity: 8dB
WS = Weight.LF(80,1,-40);   % Weight on sensitivity to assure roll-off
MU = Weight.DC(20);         % Maximum input sensitivity: 20dB
WU = Weight.HF(1e3,1,-40);  % Weight on input sensitivity to assure roll-off on the controller

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

% Design the controller. Let's start of with a simple design where the 
% weights are stacked and optimized together. This is the classical 
% way of doing H-infinity controller design
S = Channel(e/r, 'Sensitivity');
U = Channel(u/r, 'Input Sensitivity');
    
objective = [WS*S;WU*U];
constraints = [];
options = struct('controller_name','stacked\_controller');
[P,C,info] = P.solve(objective, constraints, K, options);
figure('Name','Closed loop performance'), bodemag(info,S,U)

%% 4. Discussion
% This (simple) stacked design is the classical way of doing H-infinity 
% loop shaping. The infinity norm of the combined vector is being
% minimized. However, this roughly translates to minimizing the two
% weighted transfer functions separately. This is because the 2 objectives
% are most active in different frequency regions.
%
% Why would one then switch to the 'modern' approach? We'll find out in
% example 3.
