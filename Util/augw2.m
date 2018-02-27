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

function [P] = augw2(G,W1,W2,W3)
%AUGW_IP matlab augw reimplementation for improper systems
%   AUGW_IP is the equivalent for matlabs built-in AUGW, dealing with
%   systems in the descriptor form, resulting in a proper or improper generalized
%   plant.

Wout = [];

% Setup problem
G = fromstd(G);
SysG = IOSystem(G);

systems = {SysG};

r = Signal(size(G,2));       % r = external reference to the plant (should have the dimensions of G(.,x))
u = SysG.in;           % u = controled input 
y = SysG.out;          % y = controled output
e = r - y;

connections = {};

% Add weights
if(~isempty(W1))
    W1 = fromstd(W1);
    SysW1 = IOSystem(W1);
    systems{end+1} = SysW1;
    connections{end+1} = [SysW1.in == e];
    Wout = [Wout; SysW1.out];
end
if(~isempty(W2))
    W2 = fromstd(W2);
    SysW2 = IOSystem(W2);
    systems{end+1} = SysW2;
    connections{end+1} = [SysW2.in == u];
    Wout = [Wout; SysW2.out];
end
if(~isempty(W3))
    W3 = fromstd(W3);
    SysW3 = IOSystem(W3);
    systems{end+1} = SysW3;
    connections{end+1} = [SysW3.in == y];
    Wout = [Wout; SysW3.out];
end

% solve the generalized plant
aug = IOSystem(systems{:},vertcat(connections{:}));
P = aug([Wout;e],[r;u]).content(1);
end

