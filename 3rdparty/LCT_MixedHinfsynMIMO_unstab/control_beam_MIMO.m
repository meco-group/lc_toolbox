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
addpath(genpath('C:\Users\Michiel\Dropbox\Persoonlijk\Masterthesis\cases\yalmip'))



%% Defining the subsystems
load('beam_model_MIMO')
G=beam_model_MIMO;
[A,B,C,D]=ssdata(G);
G=LTIsys(A,[B],C,D,[]);
% Weights
WS = LTIsys(Weight.LF(-60,0.5,1));
MS = LTIsys(Weight.DC(6));
WU = LTIsys(Weight.DC(10));
%% Build the plant using the LTI toolbox
lti_begin
    signal r d               % 2 new signals are created, r and u which are now available
    
    u1= G.in(1);
    u2= G.in(2);
    
    y1=G.out(3);
    y2=G.out(4);
    s1=r-G.out(1);
    s12=r-G.out(1);
    %s2=r2-G.out(2)
    K.in=[r;y1];
    K.out=[u1;u2];
%     exog_in([r1;r2])          
%     exog_out([r1-G.out(1);r2-G.out(2);u1;u2]) 
    S = Channel(s1/r,'Sensitivity');
    U1 = Channel(u1/r,'Actuator\_Effort');
     U2 = Channel(u2/r,'Actuator\_Effort2');

    
    options.gammasolver = 'mosek';
    options.orderreduction=0;
    
%     ctrl_begin('mixedHinfsyn_MIMO',options)
%         minimize(WS*S)
%         subject to
%             MS*S <= 1
%              WU*U1 <= 1
%               WU*U2 <= 1
%     ctrl_end
lti_end

Plant=extract([s1;s12;u1;u2;r;y1]/[r;u1;u2],'ol');

%% Retrieve generalized plant from LTI model
% Plant = lti_problem.gp_prob.sys_G;
my=2;
mu=2;

%%
%Control


switch 2
case 1
%Output filter
Ao=zeros(2);
Bo=[1 0 0 0 0 ;0 0 10 0 0 ];
Co=[1 0;0 1;0 0;0 0; 0 0]
Do=[0 0 0 0 0;0 0 0 0 0;0 1/10 0 0 0; 0 0 0 1/10 0;0 0 0 0 1/10];

ch.In{1}=[1 0 ];
ch.In{2}=[0 1 ];
ch.In{3}=[1 0 ];
ch.In{4}=[1 0 ];
ch.In{5}=[0 1 ];


ch.Out{1}=[1 0 0 0 0];
ch.Out{2}=[0 1 0 0 0];
ch.Out{3}=[0 0 1 0 0];
ch.Out{4}=[0 0 0 1 0];
ch.Out{5}=[0 0 0 0 1];
    case 2
        %Output filter
Ao=zeros(1);
Bo=[60 0 0 0];
Co=[1;0;0;0]
Do=[0 0 0 0;0 1/10^(4/20) 0 0; 0 0 1/10 0;0 0 0 1/10];

ch.In{1}=[1  ];
ch.In{2}=[1 ];
ch.In{3}=[1 ];
ch.In{4}=[1 ];



ch.Out{1}=[1 0 0 0];
ch.Out{2}=[0 1 0 0];
ch.Out{3}=[0 0 1 0];
ch.Out{4}=[0 0 0 1];

end


W_us=ss(Ao,Bo,Co,Do);

alpha=[1 0 0 0];

% 1.1. State-space model and dimensions
[A,B,C,D] = ssdata(balreal(Plant));
n = length(A(:,1));
mw = length(B(1,:))-mu;
mz = length(C(:,1))-my;


% 1.2. Subdivision according to inputs and outputs
                        Bw = B(:,1:mw);             Bu = B(:,mw+1:end);
Cz = C(1:mz,:);         Dzw = D(1:mz,1:mw);         Dzu = D(1:mz,mw+1:end);
Cy = C(mz+1:end,:);     Dyw = D(mz+1:end,1:mw);     Dyu = D(mz+1:end,mw+1:end);
Cout=[Cz;Cy];
Dout=[D];
beta=0;
options.gammasolver='lmilab'
options.controllersolver='lmilab'
options.beta=0;
    options.numstab=1;

options.orderreduction='on'
[K1,gamma]= mixedHinfsyn_MIMO(balreal(Plant),W_us,my,mu,alpha,ch,options)%%




%% Defining the subsystems
load('beam_model_MIMO3')
G=beam_model_MIMO3;
[A,B,C,D]=ssdata(G);
G=LTIsys(A,[B],C,D,[]);

%% Build the plant using the LTI toolbox
lti_begin
    signal r d               % 2 new signals are created, r and u which are now available
    
    u1= G.in(1);
    u2= G.in(2);
    
    y1=G.out(3);
    y2=G.out(4);
    s1=r-G.out(1);
    s12=r-G.out(1);
    %s2=r2-G.out(2)
    K.in=[r;y1];
    K.out=[u1;u2];
%     exog_in([r1;r2])          
%     exog_out([r1-G.out(1);r2-G.out(2);u1;u2]) 
    S = Channel(s1/r,'Sensitivity');
    U1 = Channel(u1/r,'Actuator\_Effort');
     U2 = Channel(u2/r,'Actuator\_Effort2');

    
    options.gammasolver = 'mosek';
    options.orderreduction=0;
    
%     ctrl_begin('mixedHinfsyn_MIMO',options)
%         minimize(WS*S)
%         subject to
%             MS*S <= 1
%              WU*U1 <= 1
%               WU*U2 <= 1
%     ctrl_end
lti_end

Plant=extract([s1;s12;u1;u2;r;y1]/[r;u1;u2],'ol');
% 1.1. State-space model and dimensions
[A,B,C,D] = ssdata(balreal(Plant));
n = length(A(:,1));
mw = length(B(1,:))-mu;
mz = length(C(:,1))-my;


% 1.2. Subdivision according to inputs and outputs
                        Bw = B(:,1:mw);             Bu = B(:,mw+1:end);
Cz = C(1:mz,:);         Dzw = D(1:mz,1:mw);         Dzu = D(1:mz,mw+1:end);
Cy = C(mz+1:end,:);     Dyw = D(mz+1:end,1:mw);     Dyu = D(mz+1:end,mw+1:end);
Cout=[Cz;Cy];
Dout=[D];
Bw=[Bw ones(length(Bw(:,1)),1)];
Dyw=[Dyw zeros(length(Dyw(:,1)),1)];
Dzw=[Dzw zeros(length(Dzw(:,1)),1)];


[Ac,Bc,Cc,Dc]=ssdata(K1);

A_cl=[A+Bu*Dc*Cy Bu*Cc;Bc*Cy Ac];
B_cl=[Bw+Bu*Dc*Dyw;Bc*Dyw];
C_cl=[Cz+Dzu*Dc*Cy Dzu*Cc];
D_cl=Dzw+Dzu*Dc*Dyw;
H=ss(A_cl,B_cl,C_cl,D_cl);
% H=series(H,W_us)

w=logspace(-5,4,100);
figure
subplot(131)
bodemag(H(1,1),w)
subplot(132)
bodemag(H(3,1),w)
subplot(133)
bodemag(H(4,1),w)

t=[0:0.001:5];
input=0.05*ones(length(t),1);
figure
y=lsim(H([1,3,4],1),input,t);
subplot(311)
plot(t,input-y(:,1))
subplot(312)
plot(t,y(:,2))
subplot(313)
plot(t,y(:,3))
figure
bode(G([1,3,4],:))