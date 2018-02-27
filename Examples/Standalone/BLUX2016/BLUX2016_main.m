clear all
close all
clc

%% 
% Example for the BENELUX MEETING ON SYSTEMS AND CONTROL 2016
% This example describes the controller design for the absolute position of
% the load of an overhead crane. It will demonstrate analysis with a
% prefixed proportional controller. This is followed by some optimal
% feedback controller design and some trade-off analysis on tracking vs
% disturbance rejection.

%% 1. System definition
g = 9.81;   %[m/s/s]
l = 0.6;    %[m]
z = 1e-3;%[-] damping of the pendulum

% G is a 2x1 lti system with as input the cart velocity and as output the
% angle and the cart position.
Gm = TFmod({[1 0];[1]},{[l 2*sqrt(g*l)*z g];[1 0]});
G = IOSystem(Gm);

% C0 is a simple proportional controller to check the performance
C0 = TFmod(100,1);
C0.name = 'proportional';
K = IOSystem(C0);

% Define some weights corresponding to what we want to achieve
WS = Weight.LF(0.4,1);
MS = Weight.DC(4);
WT = Weight.HF(5,1);
WD = Weight.DC(35);

%% 2. LTI toolbox
% Define variables to do the eventual controller design
r = Signal();    
u = G.in;

% compute the relative position from the angle
theta = G.out(1);
xr = l*theta;  

% compute the absolute position from xr and x0
x0 = G.out(2);
x = xr + x0;

% define the error on the absolute position
e = r - x;

% define the connections
conn1 = [K.in == e];               % Assign the error to the controller input
conn2 = [K.out == u];              % Assign the plant input to the controller output

% construct the closed loop system
P = IOSystem(G,K,[conn1;conn2]);

% Do the controller design
S = Channel(e/r,'Sensitivity');
U = Channel(u/r,'Input Sensitivity');
T = Channel(x/r,'Complementary Sensitivity');
D = Channel(theta/r,'Disturbance Rejection');

obj1 = WS*S;
constr1 = [MS*S <= 1, WT*T <= 1];
opts.gammasolver = 'mosek';
opts.controller_name = 'tracking';    
[P,C1,info1] = P.solve(obj1,constr1,K,opts);

obj2 = WS*S;
constr2 = [WD*D <= 1, WT*T <= 1];
opts.controller_name = 'disturbance';    
[P,C2,info2] = P.solve(obj2,constr2,K,opts);

figure, bodemag(info1,info2,S,U,T,D);
figure, step(P([x;theta;x0],r));
