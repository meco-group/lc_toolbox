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
set(cstprefs.tbxprefs,'FrequencyUnits','Hz')
addpath(genpath('..\yalmip'))
addpath(genpath('..\SDPT3-4.0'))
addpath(genpath('..\SeDuMi_1_3'))
addpath('C:\Program Files\Mosek\7\toolbox\r2013a')

addpath(genpath('C:\Users\Michiel\Documents\lti_toolbox'))

addpath(genpath('yalmip'))
addpath(genpath('SDPT3-4.0'))
addpath(genpath('../../Hinfsynontwerp/lti_toolbox'))

% Setting up the control problem
% ==============================

% Overhead crane
g = 9.81;       % gravitational acceleration
L = 0.45;       % length of pendulum
c = -0.0451;    % damping
AG = [  0    1   0 ;    % x = [theta/s, theta, xc]
       -g/L  c   0 ; 
        0    0   0 ];
BG = [0; -1/L; 1];
CG = [0, L, 1];
DG = 0;
G = ss(AG,BG,CG,DG);


%Weight values
Ms = 6;             % [dB] peak of sensitivity
WMs = tf(1, 10^(Ms/20));

Wu = 20;            % [dB] peak value of Input sensitivity
WU = tf(1, 10^(Wu/20));




%unstable output filter
Ao=0;
Bo=[1 0 0 0];
Co=[1;0;0;0];
Do=[0 0 0 0;
    0 1 0 0;
    0 0 1 0;
    0 0 0 1];
W_us=ss(Ao,Bo,Co,Do);
%augmented plant
A=AG;
Bu=BG;
Bw=[0 1;
    0 0;
    0 0]; %reference (w) and disturbance (d) as exogenous input
Cy=[0 -L -1];
Dyw=[1 0];
%position error, peak sensitivity, actuator amplification and disturbance
%rejection as design objectives
Cz=[0 -L -1;
    0 -L/10^(Ms/20) -1/10^(Ms/20);
    0 0 0;
    0 L 1];
Dzw=[1 0;
    1/10^(Ms/20) 0;
    0 0;
    0 0];
Dzu=[0;
    0;
    1/10^(Wu/20);
    0];

ch.In{1}=[1 0];
ch.In{2}=[1 0];
ch.In{3}=[1 0];
ch.In{4}=[0 1];


ch.Out{1}=[1 0 0 0];
ch.Out{2}=[0 1 0 0];
ch.Out{3}=[0 0 1 0];
ch.Out{4}=[0 0 0 1];


alpha = [1;0;0;1]; %optimize bandwidth and disturbance rejection

mu=length(Bu(1,:));
my=length(Cy(:,1));
B=[Bw Bu];
C=[Cz;Cy];
D=[Dzw Dzu;Dyw zeros(my,mu)];
P=ss(A,B,C,D);

%set solvers
options.gammasolver='lmilab' %gamma iteration
options.controllersolver='lmilab' %controller reconstruction
options.beta=0; %minimum exponential damping
options.orderreduction='on'






[K1,gamma]= mixedHinfsyn_MIMO(balreal(P),W_us,my,mu,alpha,ch,options)

%%
% Extend with state measurements to allow order reductions
Cy=[0 0 0;
    0 0 1;
    0 1 0;
    ];
    Dyw=[1 0;0 0;0 0];

my=length(Cy(:,1));
C=[Cz;Cy];
D=[Dzw Dzu;Dyw zeros(my,mu)];
P=ss(A,B,C,D);


[K2,gamma]= mixedHinfsyn_MIMO(balreal(P),W_us,my,mu,alpha,ch,options)

%%
% Some figures
close all

%build closed loops   
[Ac,Bc,Cc,Dc]=ssdata(K1);

Cy=[0 -L -1];
Dyw=[1 0];
A_cl=[A+Bu*Dc*Cy Bu*Cc;Bc*Cy Ac];
B_cl=[Bw+Bu*Dc*Dyw;Bc*Dyw];
C_cl=[Cz+Dzu*Dc*Cy Dzu*Cc];
D_cl=Dzw+Dzu*Dc*Dyw;
H01=ss(A_cl,B_cl,C_cl,D_cl);



[Ac,Bc,Cc,Dc]=ssdata(K2);

Cy=[0 0 0;
    0 0 1;
    0 1 0;
    ];
Dyw=[1 0;0 0;0 0];
A_cl=[A+Bu*Dc*Cy Bu*Cc;Bc*Cy Ac];
B_cl=[Bw+Bu*Dc*Dyw;Bc*Dyw];
C_cl=[Cz+Dzu*Dc*Cy Dzu*Cc];
D_cl=Dzw+Dzu*Dc*Dyw;
H02=ss(A_cl,B_cl,C_cl,D_cl);

figure(1);
step(H01(:,1))
hold all
step(H02(:,1))
title('step reponses')
legend('no order reduction','order reduction')


figure(2);
impulse(H01(:,2))
hold all
impulse(H02(:,2))
title('impulse disturbance response')
legend('no order reduction','order reduction')

figure(3)
bode(H01(:,1))
hold all
bode(H02(:,1))
title('reference to output')
legend('no order reduction','order reduction')

figure(4)
bode(H01(:,2))
hold all
bode(H02(:,2))
title('disturbance to outputs')
legend('no order reduction','order reduction')





%% Build the plant using the LTI toolbox
[A,B,C,D]=ssdata(G);
G=LTIsys(A,B,C,D,[]);
% Weights
WS = LTIsys(tf([1/2],[1 0]));
MS = LTIsys(Weight.DC(6));
WU = LTIsys(Weight.DC(10));
lti_begin
    signal r               % 2 new signals are created, r and u which are now available
    
    u= G.in
    y=G.out;
    s=r-G.out;
    %s2=r2-G.out(2)
    K.in=[s];
    K.out=[u];
    S = Channel(s/r,'Sensitivity');
    U1 = Channel(u/r,'Actuator\_Effort');


    
    options.gammasolver = 'lmilab';
    options.numstab=0;
    options.orderreduction=1;
    
    ctrl_begin('mixedHinfsyn_MIMO',options)
        minimize(WS*S)
        subject to
            MS*S <= 1
             WU*U1 <= 1

    ctrl_end
lti_end
