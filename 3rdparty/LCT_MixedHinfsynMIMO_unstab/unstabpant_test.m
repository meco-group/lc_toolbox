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
%%
addpath(genpath('../../lti_toolbox'))
%% build plant

Gmod = TFmod(20,conv([1 25],[1 -35]))

% Weights
WS = tf([1],[1 0]);
MS = Weight.DC(8);
WU = Weight.DC(60);


%% Solve control problem
G = IOSystem(Gmod);
K = IOSystem(1,1);
r = Signal()
u = G.in;
y = G.out;
e = r - y;

connections = [K.in == e; K.out == u];
P = IOSystem(G,K,connections);

S = Channel(e/r,'Sensitivity');
U = Channel(u/r,'Actuator\_Effort');

obj = WS*S;
constr = [MS*S<=1,WU*U<=1];
[P,C,info] = P.solve(obj,constr,K);

figure, bodemag(info)