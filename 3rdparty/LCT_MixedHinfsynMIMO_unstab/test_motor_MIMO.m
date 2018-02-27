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
% Setting up the control problem
% ==============================

load('linear_motor')
[z,p,k]=zpkdata(linear_motor,'v');
linear_motor=linear_motor+zpk(5*z,10*p,2*k)+zpk(2*z,4*p,2*k)+zpk(10*z,100*p,10*k)+zpk(z/2,p/10,k*5);
[A,B,C,D]=ssdata(linear_motor);
Bd=random('Normal',0,1,20,1);
Dd=0;
B=[B Bd];
D=[D Dd];
linear_motor=ss(A,B,C,D);
G=linear_motor;
Wt=tf([1 2*pi*10],[1 200*pi]);
Wt=Wt/freqresp(Wt,0)/2;
Wz2=tf(1,1);
Ww=tf(1,1);
%constraint
Wu=15; %dB
Wu=tf(1/10^(Wu/20),1);
Ms=6; %dB
Ms=tf(1/10^(Ms/20),1);
%Build plant
model_begin
    subsystem G Wt Wu Ms
    controller K
    signal r u y(2) e d
    
    connect to
    Gin(1) == u
    Gin(2) == d
    y(1) == r
    y(2) == Gout
    Wtin == Gout
    e == r-Gout
    e == Msin
    u == Wuin
    
    
    Kin = y
    Kout = u
    exog_in([r,d]) 
    exog_out([e;Msout;Wuout;Wtout;Gout])
model_end



Plant= lti_problem.gp_prob.sys_G;

Mz=[1 1 1 1 1]';
alpha=[1 0 0 0 1];

Ao=[0];
Bo=[2 0 0 0 0];
Co=[1;0;0;0;0];
%mz=length(Co(:,1));
Do=[0 0 0 0 0;
    0 1 0 0 0;
    0 0 1 0 0;
    0 0 0 1 0;
    0 0 0 0 1];
W_us=ss(Ao,Bo,Co,Do);

my=2;
mu=1;
Mw=[1;1];
ch.In{1}=[1 0];
ch.In{2}=[1 0];
ch.In{3}=[1 0];
ch.In{4}=[1 0];
ch.In{5}=[0 1];


ch.Out{1}=[1 0 0 0 0];
ch.Out{2}=[0 1 0 0 0];
ch.Out{3}=[0 0 1 0 0];
ch.Out{4}=[0 0 0 1 0];
ch.Out{5}=[0 0 0 0 1];

%%
options.gammasolver='mosek'
options.controllersolver='lmilab'
options.beta=0;
options.orderreduction='on'
[K,gamma]= mixedHinfsyn_MIMO6((Plant),W_us,my,mu,alpha,ch,options);%%



%%
close all


model_begin
    subsystem G 
    controller K
    signal r u y(2) e d
    
    connect to
    Gin == [u;d]
    y(1) == r
    y(2) == Gout
    e == r-Gout
        
    
    Kin = y
    Kout = u
    exog_in([r,d]) 
    exog_out([e;e;u;Gout])
model_end



Plant= lti_problem.gp_prob.sys_G;


% 1.1. State-space model and dimensions
[A,B,C,D] = ssdata(Plant);
n = length(A(:,1));
mz = 4;
mw = length(B(1,:))-mu;


% 1.2. Subdivision according to inputs and outputs
                        Bw = B(:,1:mw);             Bu = B(:,mw+1:end);
Cz = C(1:mz,:);         Dzw = D(1:mz,1:mw);         Dzu = D(1:mz,mw+1:end);
Cy = C(mz+1:end,:);     Dyw = D(mz+1:end,1:mw);     Dyu = D(mz+1:end,mw+1:end);

[Ac,Bc,Cc,Dc]=ssdata(K);

A_cl=[A+Bu*Dc*Cy Bu*Cc;Bc*Cy Ac];
B_cl=[Bw+Bu*Dc*Dyw;Bc*Dyw];
C_cl=[Cz+Dzu*Dc*Cy Dzu*Cc];
D_cl=Dzw+Dzu*Dc*Dyw;
H=ss(A_cl,B_cl,C_cl,D_cl);

%%
figure
step(H(:,1))
figure
step(H(:,2))
