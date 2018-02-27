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
clc
close all

%% 1. Connecting an augmented plant 
% Consider the system G for which we want to design a feedback controller
% K using a Hinf methodology. Two weights are being applied, MS (maximum 
% sensitivity) and WT, a robustness weight.
% Now we can start designing a controller using the lti_toolbox

%% 1.1. Let's define the systems we will be using
Gmod = ZPKmod([-2*pi*3*exp(j*pi/2.5),-2*pi*3*exp(-j*pi/2.5)],[0,-2*pi*5*exp(j*pi/2.5),-2*pi*5*exp(-j*pi/2.5)],10);

MS = Weight.DC(5); % weightDC constructs a dc weight equivalent with a peak of 5db
WS = Weight.LF(0.5,1,-60); % Construct a first order low frequency roll-off weight with a cross-over frequency of 0.5Hz
WT = zpk([-10*2*pi,-10*2*pi],[-1e3*2*pi;-1e3*2*pi],2.5e3); % robustness weight 
WU = Weight.HF(50,1,-60); % Weight on the input sensitivity to enforce roll-off in the controller

%% 1.2. Design an even better controller using the lti_toolbox
G = IOSystem(1,1);
G.add(Gmod);
K = IOSystem(1,1);  

r = Signal();
u = G.in;
y = G.out;
e = r - y;
connections = [K.in == e; K.out == u];
P = IOSystem(G,K,connections);
    
% Do the controller design
S = Channel(e/r,'Sensitivity');
U = Channel(u/r,'Input Sensitivity');
T = Channel(y/r,'Complementary Sensitivity');

objective = [];
constraints1 = [MS*S <= 1, WS*S <= 1, WT*T <= 1];
options1 = struct('controller_name','better\_controller');
[P,C1,info1] = P.solve(objective, constraints1, K, options1);

constraints2 = [MS*S <= 1, WS*S <= 1, WU*U <= 1, WT*T <= 1];
options2 = struct('controller_name','roll-off\_controller');
[P,C2,info2] = P.solve(objective, constraints2, K, options2);

figure, bodemag(info1, info2);

%% 1.3. Discussion
% This controller looks even better! Look at the high frequency part: we
% get the roll-off we desire. But what with the complementary sensitivity?
% It seems like we have some margin for improvement! Let's push the
% bandwidth to its maximum. To put it in a more technical way, let's switch
% from a feasibility problem to an optimization problem, which in fact
% means: let's minimize some function to obtain the best solution, rather
% than just picking one of the solutions! Go to level 5, ehm.. example 5 to
% see what happens.
