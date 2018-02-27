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

%% 1. Connecting an augmented plant 
% Consider the system G for which we want to design a feedback controller
% K using a Hinf methodology. Two weights are being applied, MS (maximum 
% sensitivity) and WT, a robustness weight.

%% 1.1. Let's define the systems we will be using
Gmod = ZPKmod([-2*pi*3*exp(j*pi/2.5),-2*pi*3*exp(-j*pi/2.5)],[0,-2*pi*5*exp(j*pi/2.5),-2*pi*5*exp(-j*pi/2.5)],10);

MSmod = fromstd(Weight.DC(5)); % weightDC constructs a dc weight equivalent with a peak of 5db
WTmod = ZPKmod([-10*2*pi,-10*2*pi],[-1e3*2*pi;-1e3*2*pi],2.5e3); % robustness weight 

%% 1.2. Let's now Build the generalized plant (manually) using the LCtoolbox
% Note that when you are doing a controller design, the generalized plant will be built automatically behind the scenes for you. You usually don't have to do this yourself.

% First, define a system for all systems in the plant. Add the models where you can. (note that we don't have a model for the controller K yet, so it remains empty.)
G = IOSystem(1,1);
G.add(Gmod);
K = IOSystem(1,1);
MS = IOSystem(1,1);
MS.add(MSmod);
WT = IOSystem(1,1);
WT.add(WTmod); 

% Define a bunch of signals which are interesting to you
r = Signal();               % Create some new signals, r which is our reference signal
u = G.in;
y = G.out;
e = r - y;
z = [MS.out;WT.out];

% Define which signals should be connected
connections = [MS.in == e; WT.in == y; K.in == e; K.out == u];

% Last but not least, make a new system comprised of the subsystems and the list of connections
P = IOSystem(G,K,MS,WT,connections);

%% 1.3. Retrieve the generalized plant model
% The generalized plant consists of 2 inputs; the exogeneous input r and the
% control input u, and 3 outputs; the performance outputs z and the control
% output e.
Plant = P([z;e],[r;u]).content(1);

%% 2. Building the plant via augw
% Rather than constructing the generalized plant manually, you can use
% augw_ip which internally calls the lti_toolbox. The behavior is similar
% to the augw of matlab, but allows also a connection with improper
% weights, which is beneficial for the controller design.

%% 2.1. Construct the augmented plant
Plant_augw = augw2(std(G),std(MS),[],std(WT)); %use augw2 since we use improper weights

%% 2.2. Verify the equality of both plants
if(sscomp(simplify(Plant),simplify(Plant_augw)))
    disp('Augw reference plant and augw_ip plant are equal: check passed');
else
    disp('Error! Augw reference plant and augw_ip plant are NOT equal: check NOT passed');
end
