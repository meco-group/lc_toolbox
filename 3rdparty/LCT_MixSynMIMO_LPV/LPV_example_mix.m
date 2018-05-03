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

clear all; clc; close all

k1    = SchedulingParameter('k1',[10,16],0.2);
im2 = 11;        
%im2   = SchedulingParameter('im2',[0.25,1],0.5);
param = {k1};

m1 = 10; k2 = 5; c1 = 0.2; c2 = 0.3; 
A = [0 0 1 0; 0 0 0 1; -k1*(1/m1), k1*(1/m1), -(c1+c2)/m1, c2/m1; k1*im2, -(k1+k2)*im2, c2*im2, -c2*im2];
B = [1 0 0 0; 0 1 0 0; 0 0 1/m1 0; 0 0 0 im2]; 
C = [1 0 0 0; 0 1 0 0; 0 0 0 0; 0 0 0 0; 1 0 0 0; 0 0 0 1];
D = [zeros(4,2), [0 0; 0 0; 1 0; 0 1]; zeros(2,2), zeros(2,2)];
nu = 2; ny = 2; 
sys.A = A; sys.B = B; sys.C = C; sys.D = D; sys.Ts = 0;
options = struct('var_deg',1,'var_knots',0,'controller_dependency','ada','verbose',2,'relax_deg',1);
% set up and solve the LPV control problem
        
% performance specifications

performance = 4 ;

switch performance
    case 1
        c = [0,1]; % optimize channel 2 s.t. bound on channel 1
        ch(1).performance = inf; ch(1).In = 1; ch(1).Out = [1,2];
        ch(2).performance = inf; ch(2).In = 2; ch(2).Out = [3,4];
    case 2
        c = [0,1]; % optimize channel 2 s.t. bound on channel 1
        ch(1).performance = inf; ch(1).In = 1; ch(1).Out = [1,2];
        ch(2).performance = 2; ch(2).In = 2; ch(2).Out = [3,4];
    case 3
        c = [1,1]; % bounds on both channels
        ch(1).performance = 2; ch(1).In = 1; ch(1).Out = [1,2];
        ch(2).performance = 2; ch(2).In = 2; ch(2).Out = [3,4];
    case 4
        c = [1,1]; % bounds on both channels
        ch(1).performance = inf; ch(1).In = 1; ch(1).Out = [1,2];
        ch(2).performance = inf; ch(2).In = 2; ch(2).Out = [3,4];
    otherwise
        error('Wrong selection');
end

[K,sol_info] = LPV_unstructured_mix(sys,ny,nu,c,ch,param,options)