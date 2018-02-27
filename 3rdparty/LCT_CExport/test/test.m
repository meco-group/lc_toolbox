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

settings.namespace = 'Ballbot';
settings.filename = 'attitude_controller';
settings.cname = 'AttitudeController';

Ts = 0.01;

A1 = [-0.1916, -0.1814; 0.1814, -0.1916];
B1 = [0; 0.5447];
C1 = [-0.1379, 0.6199];
D1 = [0];
c1 = ss(A1,B1,C1,D1,Ts);

A2 = [-0.3573, 0.2659, -0.3172; 0.2523, -0.09952, -0.5243; -0.3281, -0.5175, 0.131];
B2 = [ 0; 0.06922; 2.487];
C2 = [-1.666, -0.4159, -0.08422];
D2 = [0];

c2 = ss(A2,B2,C2,D2,Ts);

n1 = 'controller1';
n2 = 'controller2';

c_export(settings,c1,n1,c2,n2);