clear all
close all
clc

%% Controller design for a ball balancing robot - attitude loop
% The zero in zero in the attitude loop is cancelled by the PI controller
% which acts as a prefilter to the designed controller.

% Declare model
Gm = ZPKmod(0,[-0.9,0.9],4e-3);
PIm = ZPKmod(-5,0,1);

G = IOSystem(Gm);
PI = IOSystem(PIm);
K = IOSystem(1,1);

th_ref = Signal();
u = G.in;
th = G.out;
e = th_ref - th;

connections = [K.in == e; K.out == PI.in; PI.out == u];
P = IOSystem(G,PI,K,connections);    

% Declare shape functions ad channels
WS = Weight.LF(2,1,-10);
MS = Weight.DC(6);
WU = Weight.HF(10,1,-30)*1e-3;
    
S = Channel(e/th_ref,'Sensitivity');
U = Channel(u/th_ref,'Input Sensitivity');
T = Channel(th/th_ref,'Compl. Sensitivity')

% Specify the optimal controller design problem
obj = WU*U;
constr = [MS*S <= 1, WS*S <= 1];
[P,C,info] = P.solve(obj,constr,K);
figure, bodemag(info,S,U,T)

% Clean up the controller and plot the results
C = removehigh(ssminreal(C),1e3);
C.name = 'clean';
K.add(C)
figure, bode(K.content(1)*PIm,K.content(2)*PIm)
