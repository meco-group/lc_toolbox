clear all
close all
clc

addpath(genpath('../lti_toolbox'));

%% Satelite Tracking example
% The system describes a controller design for a disc-antenna which is to
% track a satelite oribiting through space. Zero steady-state error under 
% constant is disturbances is desirable.
% Furthermore, a bandwidth of 0.1Hz is needed. The remaining freedom should
% be used to minimize the actuation effort, which also makes the controller
% robust to high-frequency noise. In that regard, some roll-off on the
% controller is also desirable

%% 1.1. Define the system and select some weights
mod1 = TFmod(10,[1 15 50 0]);
%mod2 = FRDmod(mod1,logspace(1e-1,1e2,50),'Hz');

G = IOSystem(1,1);
G.add(mod1);
K = IOSystem(1,1);
r = Signal();

e = r - G.out;
u = G.in;
y = G.out;
conn = [K.in == e; K.out == G.in];
CL = IOSystem(G,K,conn);

WS = Weight.LF(0.1,2,-40);
MS = Weight.DC(4);
MS2 = Weight.DC(8);
WU = Weight.HF(3,1,-40);
WT = Weight.HF(1,3);

%% 1.2. Design the optimal controller using the lti_toolbox
S = Channel(e/r,'Sensitivity');
U = Channel(u/r,'Input sens.');
T = Channel(y/r,'Compl. sens.');

[CL,C1,info1] = CL.solve(WT*T,[WS*S <= 1,MS*S <= 1],K);
[CL,C2,info2] = CL.solve(WT*T,[WS*S <= 1,MS2*S <= 1],K);

figure, bodemag(info1,info2,S,T)

%% 1.3. Time domain simulation
u = @(t) (t>0);
t = 20;
sim(T(CL),u,t);

%%
[CL,C3,info3] = CL.solve(WU*U,[WS*S <= 1,MS*S <= 1],K);