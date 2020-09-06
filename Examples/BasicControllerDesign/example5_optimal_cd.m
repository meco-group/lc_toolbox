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

%% 1. Connecting an augmented plant 
% Consider the system G for which we want to design a feedback controller
% K using a Hinf methodology. Two weights are being applied, MS (maximum 
% sensitivity) and WT, a robustness weight.
% Now we can start designing a controller using the lti_toolbox

%% 1.1. Let's define the systems we will be using
Gmod = ZPKmod([-2*pi*3*exp(j*pi/2.5),-2*pi*3*exp(-j*pi/2.5)],[0,-2*pi*5*exp(j*pi/2.5),-2*pi*5*exp(-j*pi/2.5)],10);

MS = Weight.DC(5); % weightDC constructs a dc weight equivalent with a peak of 5db
WS = Weight.LF(0.01,2,-40); % Construct a first order low frequency roll-off weight with a cross-over frequency of 0.5Hz
WT = zpk([-10*2*pi,-10*2*pi],[-1e3*2*pi;-1e3*2*pi],2.5e3); % robustness weight 
WU = Weight.HF(50,1,-60); % Weight on the input sensitivity to enforce roll-off in the controller

WSu = Weight.LF(0.01,2); % Unstable variant of WS
WUi = Weight.HF(50,1); % Improper variant of WU

%% 1.2. Design the optimal controller using the lti_toolbox
% Let's also set some options to get the plots we want
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

% objective1 = [];
% constraints1 = [MS*S <= 1, WS*S <= 2, WU*U <= 1, WT*T <= 1];
% options1 = struct('controller_name','roll-off\_controller');
% [P,C1,info1] = P.solve(objective1, constraints1, K, options1);

objective2 = WS*S;
constraints2 = [MS*S <= 1; WU*U <= 1];
options2 = struct('controller_name','optimal\_controller');
[P,C2,info2] = P.solve(objective2, constraints2, K, options2);

objective3 = WSu*S;
constraints3 = [MS*S <= 1; WUi*U <= 1];
options3 = struct('controller_name','optimal\_controller');
[P,C3,info3] = P.solve(objective3, constraints3, K, options3);

% figure, showall(info1, info2, info3);
figure, showall(info2,info3)

%% 1.3. Discussion
% I think it is quite clear that the performance has increased! Some things
% to look at:
% 1: The bandwidth (sensitivity!) has increased which means we have better
% disturbance rejection.
% 2: This can also be seen in the controller: the (dc) gain of the
% controller has increased: 10dB to 20dB
% 3: Which constraint is now active? Looking at S,U and T, it becomes clear
% that the constraint on U is now limiting the attainable bandwidth: better
% tracking requires more effort of the controller, but we limited this,
% remember? 
%
% Improving the design: 
% 1: If the constraint on U was used to enforce roll-off, but without fair
% knowledge of the actual actuator constraints, the constraint can be put a
% relaxed as to allow an even larger bandwidth. The best you can attain of
% course is to get near WT, since this ensures robustness of the
% controller!
% 2: Cleaning up the controller will result in a flat controller near low
% frequencies. However, typically we want something like an integrator in
% our controller. Why is it not in there at the moment? Because there is a
% controller in the plant! This integrator ensures already -20dB roll-off
% on the sensitivity function. So if we impose a -40dB roll-off weight, we
% will attain an integrator in the controller!