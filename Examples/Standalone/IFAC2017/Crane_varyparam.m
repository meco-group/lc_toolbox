clear all
close all
clc

%% 
% Example for the IFAC world congress 2017
% This example describes the controller design for the absolute position of
% the load of an overhead crane. The model is depending on the pendulum's
% length l. Small angles are assumed so that we obtain an LPV system.

%% 1. System definition
% 1.1. some fixed parameters
g = 9.81;   %[m/s/s]
z = 1e-3;   %[-] damping of the pendulum

% 1.2. varying parameter
% The model we are implementing has the cart's speed as an input and the
% absolute load position on the output. The system is as follows:
% X/U = tf([2s^2 ds/l g/l],[s^3 ds^2/l gs/l])
% Since the parameter appears a fraction in the according state space, the
% parameter 'a' is chosen as the inverse of l, so that:
% X/U = tf([2s^2 das ga],[s^3 das^2 gas])
% The bounds on the new parameter can be computed based on the bounds for l
Lbound = [0.2,0.8]; Lrate = [-0.1,0.1];

% 1.3. The system parameter
l = SchedulingParameter('L',Lbound,Lrate);    %[m]
lg = linspace(Lbound(1),Lbound(2),100);

% Linear approximation
degree = 1; knots = 2;
opti = splines.OptiSpline();
b = splines.BSplineBasis(Lbound,degree,2+knots);
t = splines.TensorBasis(b,'L');
linv = opti.Function(t,[1,1],'full');

grid = linspace(Lbound(1),Lbound(2),100)';
var = linv.list_eval(grid);
ref = 1./grid;

opti.minimize(norm(ref-var,2));
opti.subject_to({var-ref >= 0});
opti.solver('ipopt');
sol = opti.solve();
linv = sol.value(linv); %inverse mapping of l
figure, plot(grid,linv.list_eval(grid),grid,ref);

% 1.4 The model and the system
A = [0 1 0;0 0 1;0 -g*linv -z*linv];
B = [0;0;1];
C = [g*linv,z*linv,2];
D = [0];

Gm = SSmod(A,B,C,D,l);
G = IOSystem(Gm)

% Define some weights corresponding to what we want to achieve
WS = Weight.LF(0.10,1,-20);
MS = Weight.DC(6);
WU = Weight.HF(0.1,1,-40);

%% 2. LTI toolbox
% controller system
K = IOSystem(1,1);

% Define variables to do the eventual controller design
r = Signal();

% assign the convenience variable u
u = G.in;
x = G.out;

% define the error on the absolute position
e = r - x;

c1 = (K.in == e);               % Assign the error to the controller input
c2 = (K.out == u);              % Assign the plant input to the controller output

P = IOSystem(G,K,[c1;c2]);

% Do the controller design
S = Channel(e/r,'Sensitivity');
U = Channel(u/r,'Input Sensitivity');
T = Channel(x/r,'Complementary Sensitivity');

options.controller_name = 'Varying parameter controller';
obj = []; %[WS*S;WU*U];
constr = [WS*S];
[P,C,info] = P.solve(obj,constr,K,options);
figure, bodemag(info);

%% 3. Time domain simulation
simulation_on = true;
if simulation_on
    P = extract(T);

    t = [0,150];
    x0 = zeros(P.nx,1);
    
    % Plot the parameter domain
    A = 0.1;
    f = 1/(2*pi);
    p = @(t) (0.5+A*sin(f*2*pi*t)); %fixed parameter
    dp = @(t) (A*cos(f*2*pi*t)*2*pi*f);
    tp = linspace(0,150,1500);
    
    figure, 
    x = [Lbound(1),Lbound(1),Lbound(2),Lbound(2)];
    y = [Lrate(1),Lrate(2),Lrate(2),Lrate(1)];
    patch(x,y,'blue','FaceAlpha',.3), hold on
    plot(p(tp),dp(tp),'k','LineWidth',2), hold off
    axis([Lbound(1)-0.05,Lbound(2)+0.05,Lrate(1)-0.025,Lrate(2)+0.025])
    xlabel('$L$','interpreter','latex')
    ylabel('$\dot{L}$','interpreter','latex')

    u = @(t) (0.5*sin(0.05*2*pi*t)); %step input
    [y,t,x] = sim(P(x,r),u,p,t);
    
    
%     % Do the simulation
%     opts = odeset('AbsTol',1e-3);
%     [t,x] = ode23(@myODE,t,x0,opts,P,p,u);
% 
%     N = length(t);
%     y = zeros(1,N);
%     for k = 1:N
%         y(1,k) = myOutput(t(k),x(k,:).',P,p,u);
%     end
    figure, plot(t,x)
    figure, plot(t,y,'b',t,u(t),'k')

%     save('crane_varyparam_sim')
end