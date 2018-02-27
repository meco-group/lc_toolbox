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

classdef ControllerDesign < ControlProblem
    %CONTROLLERDESIGN Supports controller design
    %   Controller design class. Stores a set of constraints and objectives
    %   and calculates the optimal solution when executed.
    %   Performances can easily be checked
    
    properties
        method = 'none';    % Solver to use: mixedHinfsyn, compute_FO_mix, compute_RO_mix, none
        order = -1;         % Don't include order constraint to begin with
        
        solver = [];             % solver object
        processor = SimpleProcessor();
    end
    
    methods
        function obj = ControllerDesign(varargin)
            obj@ControlProblem(varargin{:});
            obj.options = struct('FullOrderSolver',[],'ReducedOrderSolver',[]);
        end
    end
end

