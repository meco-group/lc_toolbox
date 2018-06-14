%% This is an example to test the discrete time LPV solver. 
% Taranjitsingh Singh : 17.11.2017

% This example is taken from the paper : Masayuki Sato : "Discrete-time
% Gain-Scheduled Output-Feedback Controllers Exploiting Inexact Scheduling
% Parameters via Parameter-Dependent Lyapunov Functions"

% This scripts validates for LPV discrete time examples for various cases. 

 Ts = 1/100; %100 Hz sampling frequency
 mu = 0.25;
  
 example_type = 1;
 switch example_type
     case 1
          p = SchedulingParameter('p',[0,1],0);
          A = mu*[1-p    0           -2+p;
                 2-p     -1          1-p;
                 -1+p    1-(3*p)     -p];
          param = {p};
     case 2
         p = SchedulingParameter('p',[0,1],0.2);
           A = mu*[1-p    0           -2+p;
         2-p     -1          1-p;
         -1+p    1-(3*p)     -p];
          param = {p};
     case 3
         p = SchedulingParameter('p',[0,1],0);
         t = SchedulingParameter('t',[0,1],0);
            A = mu*[1-p+t    0           -2+p;
         2-p     -1          1-p;
         -1+p    1-(3*p)     -p];
          param = {p,t};
     case 4
         p = SchedulingParameter('p',[0,1],0.1);
         t = SchedulingParameter('t',[0,1],0);
            A = mu*[1-p+t    0           -2+p;
         2-p     -1          1-p;
         -1+p    1-(3*p)     -p];
          param = {p,t};
     case 5
         p = SchedulingParameter('p',[0,1],0.1);
         t = SchedulingParameter('t',[0,1],0.2);
            A = mu*[1-p+t    0           -2+p;
         2-p     -1          1-p;
         -1+p    1-(3*p)     -p];
          param = {p,t};         
     otherwise
         error('Wrong selection of example type');
 end

 B = [1;0;0];
 
 C = [1 0 0];
 
 D = 0;
 
 G = DSSmod(A,B,C,D,eye(3),param,Ts);

H = IOSystem(1,1);
H.add(G);

K = IOSystem(1,1);

WS = c2d(Weight.LF(0.4,1,-40),Ts); % for sensitivity
WU = 0.8*c2d(Weight.HF(1,1,-5),Ts); % for control sensitivity
MS = c2d(Weight.DC(2),Ts);

r = Signal();u = H.in;y = H.out;e = r - y;

connections = [u==K.out;e==K.in];
CL = IOSystem(H,K,connections);

S = Channel(e/r,'Sensitivity');
U = Channel(u/r,'Actuator Effort');
T = Channel(y/r,'Compl Sensitivity');

synthesis_type = 2;

switch synthesis_type 
    case 1
        obj = [WS*S];
        constr = [];
    case 2
        obj = [WS*S];
        constr = [WU*U<=1, MS*S <= 1];
    otherwise
        error('Wrong selection of synthesis type');
end

options = struct('var_deg',1,'var_knots',0,'verbose',2,'controller_dependency','a');
[CL,K,sol_info] = CL.solve(obj,constr,K,options);

% S = CL(e,r);
% U = CL(u,r);
% gam = sol_info.gamma(1);
% figure; bodemag([S;U],[gam*WS^-1;WU^-1]);
mod = CL.model();
subplot(1,3,1)
sigma(mod(e,r))
hold on
sigma(sqrt(sol_info.gamma(1))*inv(WS))
sigma(inv(MS),'k');
title ('Sensitivity')
subplot(1,3,2)
sigma(mod(u,r))
hold on
sigma(inv(WU))
title ('Control Sensitivity')
subplot(1,3,3)
step(mod(y,r))
title ('Step Response')