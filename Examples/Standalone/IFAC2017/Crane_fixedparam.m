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
z = 1e-4;   %[-] damping of the pendulum

% 1.2. varying parameter
% The model we are implementing has the cart's speed as an input and the
% absolute load position on the output. The system is as follows:
% X/U = tf([2s^2 ds/l g/l],[s^3 ds^2/l gs/l])
% Since the parameter appears a fraction in the according state space, the
% parameter 'a' is chosen as the inverse of l, so that:
% X/U = tf([2s^2 das ga],[s^3 das^2 gas])
% The bounds on the new parameter can be computed based on the bounds for l
Lbound = [0.2,0.8]; Lrate = [0,0];
% Abound = 1/Lbound since l>0:
Abound = sort(1./Lbound)
% Arate = -(1/l^2)*dl/dt. So the rate of a lies within the interval:
Arate = sort(-max(Abound.^2)*Lrate)

% 1.3. The system parameter
a = SystemParameter(Abound,Arate); %[1/m]

% 1.4 The system
A = [0 1 0;0 0 1;0 -g*a -z*a];
B = [0;0;1];
C = [g*a,z*a,2];
D = [0];
G = SSsys(A,B,C,D);

% Define some weights corresponding to what we want to achieve
WS = Weight.LF(0.1,2,-40);
WU = Weight.HF(0.1,1,-40);

%% 2. LTI toolbox
options.output = struct('controller',1,'performance',1,'closedloop',0,'displaystyle','conference');

lti_begin(options)
    % Define variables to do the eventual controller design
    r = Signal();
    
    % assign the convenience variable u
    u = G.in;
    x = G.out;
    
    % define the error on the absolute position
    e = r - x;

    K.in = e;               % Assign the error to the controller input
    K.out = u;              % Assign the plant input to the controller output

    % Do the controller design
    S = Channel(e/r,'Sensitivity');
    U = Channel(u/r,'Input Sensitivity');
    T = Channel(x/r,'Complementary Sensitivity');
    show(S,T,U)
    
    ctrl_begin('fixed parameter controller');
        minimize([WS*S;WU*U])
    ctrl_end
lti_end

%% 2. Time domain simulation
simulation_on = false;
if simulation_on
    P = extract(T);

    t = [0,200];
    x0 = zeros(P.nx,1);

    A = 0.1;
    f = 1/(2*pi);
    p = @(t) (1/(0.5+A*sin(f*2*pi*t))); %fixed parameter
    dp = @(t) (A*cos(f*2*pi*t)*2*pi*f); 

    u = @(t) (0.5*sin(0.05*2*pi*t)); %step input
    [t,x] = ode45(@myODE,t,x0,[],P,p,u);
    N = length(t);
    y = zeros(1,N);
    for k = 1:N
        y(1,k) = myOutput(t(k),x(k,:).',P,p,u);
    end
    figure, plot(t,x)
    figure, plot(t,y,'b',t,u(t),'k')
end

%% 3. Solve local problems
solve_lti = false;
if solve_lti
    Ge = G.f(5); % evaluate 5 local equidistant problems

    for k = 1:length(Ge)
        G = Ge{k};

        clear lti_problem
        lti_begin(options)
            % Define variables to do the eventual controller design
            r = Signal();

            % assign the convenience variable u
            u = G.in;
            x = G.out;

            % define the error on the absolute position
            e = r - x;

            K.in = e;               % Assign the error to the controller input
            K.out = u;              % Assign the plant input to the controller output

            % Do the controller design
            S = Channel(e/r,'Sensitivity');
            U = Channel(u/r,'Input Sensitivity');
            T = Channel(x/r,'Complementary Sensitivity');
            show(S,U,T)

            ctrl_begin('tracking');
                minimize([WS*S;WU*U])
            ctrl_end
        lti_end
    end
end
