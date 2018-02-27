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
clear global
clc
set(cstprefs.tbxprefs,'FrequencyUnits','Hz')
addpath(genpath('C:\Users\Michiel\Dropbox\Persoonlijk\Masterthesis\cases\yalmip'))
addpath(genpath('C:\Users\Michiel\Dropbox\Persoonlijk\Masterthesis\cases\SDPT3-4.0'))
addpath(genpath('C:\Users\Michiel\Dropbox\Persoonlijk\Masterthesis\cases\SeDuMi_1_3'))
addpath('C:\Program Files\Mosek\7\toolbox\r2013a')
addpath(genpath('C:\Users\Michiel\Dropbox\Persoonlijk\Masterthesis\cases\lti_toolbox'))
%% Defining the subsystems
load('beam_model')
G=beam_model;



%% Build the plant using the LTI toolbox
model_begin
    subsystem G     % G, 
    controller K            % K is the definition of the controller
    signal r1 u1       % 2 new signals are created, r and u which are now available
    
    connect to
    y = Gout(2);
    u1 == Gin                % u is connected to the first input of G
    Kin = [r1;y]                 % e is set as the control output for the generalized plant (Can also be inputed as 'control_out([e])')
    Kout = [u1]                % u is set as the control input for the generalized plant (Can also be inputed as 'control_in([u])')
    exog_in([r1])            % r is an exogeneous input to the system 
    exog_out([r1-Gout(1);r1-Gout(1);u1])     % W1.y(1) and W2.y(1) are exognenous outputs for the system
model_end

%% Retrieve generalized plant from LTI model
Plant = lti_problem.gp_prob.sys_G;
my=2;
mu=1;
%Control

%constraint
W_u=10 %dB
W_s=6 %dB
alpha=[1 0 0 ];

ch.In{1}=[1];
ch.In{2}=[1];
ch.In{3}=[1];


ch.Out{1}=[1 0 0];
ch.Out{2}=[0 1 0];
ch.Out{3}=[0 0 1];

%Output filter
Ao=0;
Bo=[100 0 0];
Co=[1;0;0];
Do=[0 0 0;0 0 1/10^(W_u/20) ;0 1/10^(W_s/20) 0 ];
W_us=ss(Ao,Bo,Co,Do);

% 1.1. State-space model and dimensions
[A,B,C,D] = ssdata(Plant);
n = length(A(:,1));
mw = length(B(1,:))-mu;
mz = length(C(:,1))-my;


% 1.2. Subdivision according to inputs and outputs
                        Bw = B(:,1:mw);             Bu = B(:,mw+1:end);
Cz = C(1:mz,:);         Dzw = D(1:mz,1:mw);         Dzu = D(1:mz,mw+1:end);
Cy = C(mz+1:end,:);     Dyw = D(mz+1:end,1:mw);     Dyu = D(mz+1:end,mw+1:end);
Cout=[Cz;Cy];
Dout=[D];

%%
options.gammasolver='lmilab'
options.controllersolver='lmilab'
options.beta=0;
options.orderreduction='on'
[K1,gamma]= mixedHinfsyn_MIMO(balreal(Plant),W_us,my,mu,alpha,ch,options);%%


%%
[Ac,Bc,Cc,Dc]=ssdata(K1);

A_cl=[A+Bu*Dc*Cy Bu*Cc;Bc*Cy Ac];
B_cl=[Bw+Bu*Dc*Dyw;Bc*Dyw];
C_cl=[Cz+Dzu*Dc*Cy Dzu*Cc];
D_cl=Dzw+Dzu*Dc*Dyw;
H=ss(A_cl,B_cl,C_cl,D_cl);

step(H)

figure
bode(H)

