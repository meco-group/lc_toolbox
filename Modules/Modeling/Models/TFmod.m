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

function sys = TFmod(varargin)
% Creates a model based on transfer function data.
% LCToolbox counterpart of MATLAB's \c tf().


if (all(cellfun(@isnumeric,varargin)) || any(cellfun(@iscell,varargin)))
tfsys = tf(varargin{:});
[A,B,C,D,E,Ts] = dssdata(tfsys);
sys = DSSmod(A,B,C,D,E,Ts);
else
    if nargin <4 
        Ts = 0;
    elseif nargin < 3
        error('Not enough arguments');
    else
        Ts = varargin{4};
    end
    num = varargin{1};
    den = varargin{2};
    param = varargin{3};
    sys = tf2sspar(num,den,param,Ts);
end

% % TODO: spline implementation
% switch nargin
%     case 1
%         % static gain
%         
%     case 2
%         % ct zpk
%         
%     case 3
%         % dt zpk
%         
% end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [G] = tf2sspar(num,den,param,Ts)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

[mnum,nnum] = size(num);
[mden,n] = size(den);
A = [];B = [];C = []; D = [];
param = param;
Ts = Ts;

if mnum > 1
    num = num';
end
if mden > 1
    den = den';
end
if isnumeric(den(1))
num = [zeros(mnum,n-nnum) num]/den(1);
else
num = [zeros(mnum,n-nnum) num]/den(1).coeff.data(1);
end
if ~isempty(num)
    D = num(:,1);
else
    D = [];
end
if isnumeric(den(1))
den = den(2:n)/den(1);
else
den = den(2:n)/den(1).coeff.data(1);
end
A = [-den; eye(n-2,n-1)];
B = eye(n-1,1);
if mnum > 0
    C= num(:,2:n) - num(:,1) * den;
else
    C = [];
end
G = DSSmod(A,B,C,D,eye(size(A,1)),param,Ts);
end