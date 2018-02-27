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

function mod = SSmod(varargin)
% Creates a model based on state-space data.
% LCToolbox counterpart of MATLAB's \c ss().

if isa(varargin{1},'Gridmod')
    [~,~,~,~,~,mod] = smile_techniqueOC_spline(varargin{:});
else
    if nargin >= 4
        A = varargin{1};
        B = varargin{2};
        C = varargin{3};
        D = varargin{4};
        varargin(1:4) = [];    
    elseif nargin >= 1
        D = varargin{1};
        A = zeros(0,0);
        B = zeros(0,size(D,2));
        C = zeros(size(D,1),0);
        varargin(1) = [];
    end

    E = eye(size(A));
    mod = DSSmod(A,B,C,D,E,varargin{:});
end
end

