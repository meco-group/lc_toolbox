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

function self = IOSystem(varargin)
% Creates a SystemOfModels or SystemOfSystems object, depending on the
% input arguments. The end user is not supposed to construct such objects
% himself but only through this function. 
% 
% Parameters:
%  varargin: may contain a mix of different Model objects, SystemOfModels
%  objects or SystemOfSystems objects. The last input argument, furthermore, 
%  is allowed to be a cell that specifies how all models and/or systems are 
%  connected to each other. 
%
% Return values:
%  self : a SystemOfModels or a SystemOfSystems object @type AbstractSystem

if nargin > 0
    if isa(varargin{1},'AbstractSystem')
        self = SystemOfSystems(varargin{:});
    else
        self = SystemOfModels(varargin{:});
    end
end