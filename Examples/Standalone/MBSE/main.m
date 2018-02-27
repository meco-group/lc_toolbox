clear all; close all; clc;

data = load('sys.mat');

Gmodel = fromstd(data.Sd);
Ts = Gmodel.Ts;

G = IOSystem(d2c(Gmodel));      % Make a system for the plant
K = IOSystem(2,2);      % Make a system for the controller

% Weight to maximize the bandwidth
%WS = c2d(Weight.LF(0.5,1,-40),Ts,'Tustin'); % sensitivity weight

W1 = Weight.DC(5);

WS = [Weight.DC(5) Weight.DC(5);Weight.DC(5) Weight.DC(5)];
WT = [Weight.DC(5) Weight.DC(5);Weight.DC(5) Weight.DC(5)];

% Make the control configuration
r = Signal(2); % reference input
u = G.in;     % control input
y = G.out;    % measured output
e = r - y;    % error on the absolute position

connections = [K.in == e;K.out == u];
CL = IOSystem(G,K,connections);

% Controller design
S1 = Channel(e(1)/r(1),'Sensitivity');
T = Channel(y/r,'Compl. Sensitivity');

% Controller design 1
obj = [W1*S1] %WS*S];
constr = [] %WT*T <= 1];
CL.solve(obj,constr,K);

% Compute closed loop response
mod = CL.model();
