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

function [num, den] = butter2(n, Wn, varargin)
%BUTTER2 Compute continuous or discrete time butterworth filter
%   n = order
%   Wn = cut-off frequency [rad/s]
%   Ts (optional) = sample time
%   type = 'low' or 'high' ('stop and pass are undefined')

type = 'low';
analog = false;

for k = 1:length(varargin)
    if length(varargin{k}) == 1
        analog = strcmp(varargin{k},'s');
    else
        type = varargin{k};
    end
end
    
if ~analog,
    fs = 2;
    u = 2*fs*tan(pi*Wn/fs);
else
    u = Wn;
end

% continuous time butterworth filter
if strcmp(type,'high')
    num = [1 zeros(1,n)]; 
else
    num = u^n;   
end

den = butterpoly(n,u);

% Discretize if necessary
if ~analog
    [num,den] = tfdata(c2d(tf(num,den),1/fs,'Tustin'),'v');
end
end