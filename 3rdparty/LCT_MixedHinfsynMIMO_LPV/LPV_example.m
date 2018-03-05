
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

% load stuff
meco_binaries('cpp_splines','develop');
import splines.*

example_type = 'uvar_bounded';

%scheduling parameters
switch example_type
    case 'uvar'
        k1    = SchedulingParameter('k1',[10,16],0);
        im2   = 0.625;
        param = {k1}; % stack SchedulingParameters in a cell for solver
    case 'uvar_bounded'
        k1    = SchedulingParameter('k1',[10,16],[-2,2]);
        im2   = 0.625;
        param = {k1}; % stack SchedulingParameters in a cell for solver
    case 'uvar_unbounded'
        k1    = SchedulingParameter('k1',[10,16],[-inf,inf]);
        im2   = 0.625;
        param = {k1}; % stack SchedulingParameters in a cell for solver
    case 'mvar'
        k1    = SchedulingParameter('k1',[10,16],0);
        im2   = SchedulingParameter('im2',[0.25,1],0);
        param = {k1,im2};
    case 'mvar_k1_bounded'
        k1    = SchedulingParameter('k1',[10,16],[-2,2]);
        im2   = SchedulingParameter('im2',[0.25,1],0);
        param = {k1,im2};  
    case 'mvar_im2_bounded'
        k1    = SchedulingParameter('k1',[10,16],0);
        im2   = SchedulingParameter('im2',[0.25,1],[-0.05,0.05]);
        param = {k1,im2};
    otherwise
        error('Wrong selection of example_type . . .');
end

% system [Interpolation-Based Modeling of MIMO LPV Systems, J. De Caigny, J. F. Camino, and J. Swevers, IEEE TCST 2011]
m1 = 10; k2 = 5; c1 = 0.2; c2 = 0.3; 
A = [0 0 1 0; 0 0 0 1; -k1*(1/m1), k1*(1/m1), -(c1+c2)/m1, c2/m1; k1*im2, -(k1+k2)*im2, c2*im2, -c2*im2];
B = [1 0 0 0; 0 1 0 0; 0 0 1/m1 0; 0 0 0 im2]; 
C = [1 0 0 0; 0 1 0 0; 0 0 0 0; 0 0 0 0; 1 0 0 0; 0 0 0 1];
D = [zeros(4,2), [0 0; 0 0; 1 0; 0 1]; zeros(2,2), zeros(2,2)];
nu = 2; ny = 2; 
sys.A = A; sys.B = B; sys.C = C; sys.D = D; sys.Ts = 0;
options = struct('var_deg',1,'var_knots',0,'controller_dependency','a','spec',inf,'objective','wc','verbose',2);
% set up and solve the LPV control problem
[K,sol_info] = LPV_unstructured_OF(sys,ny,nu,param,options)

%%
k1 = 13;
im2 = 0.625;
A = [0 0 1 0; 0 0 0 1; -k1*(1/m1), k1*(1/m1), -(c1+c2)/m1, c2/m1; k1*im2, -(k1+k2)*im2, c2*im2, -c2*im2];
B = [1 0 0 0; 0 1 0 0; 0 0 1/m1 0; 0 0 0 im2]; 
C = [1 0 0 0; 0 1 0 0; 0 0 0 0; 0 0 0 0; 1 0 0 0; 0 0 0 1];
D = [zeros(4,2), [0 0; 0 0; 1 0; 0 1]; zeros(2,2), zeros(2,2)];
sys.A = A; sys.B = B; sys.C = C; sys.D = D; sys.Ts = 0.0001;
gen_sys = extract_generalized_plant(sys,nu,ny);
A  = gen_sys.A;  Bw  = gen_sys.Bw;  Bu  = gen_sys.Bu;
Cz = gen_sys.Cz; Dzw = gen_sys.Dzw; Dzu = gen_sys.Dzu;
Cy = gen_sys.Cy; Dyw = gen_sys.Dyw; Dyu = zeros(nu,ny);
nx = gen_sys.nx; nz = gen_sys.nz; nw = gen_sys.nw;
opts.Ts = 0.0001;
opts.Dc = 0;
Info = compute_FO_Hinf_LTI(A,Bw,Bu,Cz,Cy,Dzw,Dzu,Dyw,Dyu,opts)


