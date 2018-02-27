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

tfsys = tf(varargin{:});
[A,B,C,D,E,Ts] = dssdata(tfsys);
sys = DSSmod(A,B,C,D,E,Ts);

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

