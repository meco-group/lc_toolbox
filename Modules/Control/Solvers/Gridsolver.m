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

classdef Gridsolver < Solver
    % Gridsolver allows to easily solve a grid of control problems, for
    % example to compare a gain-scheduled controller with a parametrized
    % LPV controller. 
    
    properties
        cell_solver     % cell of appropriate Solver objects @type cell
    end
    
    methods
        function self = Gridsolver(cell_solver,options)
        % Constructor for Gridsolver objects.
        % 
        % Parameters:
        %  cell_solver: the grid of Solver objects @type cell
        %  options: options that the user wants to pass to the solvers of the grid elements @type struct
        %
        % Return values:
        %  self: the Gridsolver object @type Gridsolver
        
            self.cell_solver = cell_solver;
            self.options = options;
        end
        
        function self = solve(self,config,specs,vars)
        % Parses the control problem for every subproblem in the grid such
        % that it can be solved.
        % 
        % Parameters:
        %  self: the solver interface @type Gridsolver
        %  config: the control configuration @type SystemOfSystems
        %  specs: specifications of the control problem @type
        %  ControllerDesign
        %  vars: \c cell containing the optimization variables @type cell
        %
        % Return values:
        %  self: the Solver object (hopefully) containing the solutions for every subproblem @type Gridsolver
        
            function [K,H,gamma,mu,performance,solved,info] = single_solve(solver,model,specs,vars)
                solution = solve(solver,model,specs,vars);
                K = solution.K;
                H = solution.H;
                gamma = solution.gamma;
                mu = solution.mu;
                performance = solution.performance;
                solved = solution.solved;
                info = solution.info;
            end
            
            cconf = content(remove_empty(unpack(config)));
            isg = cellfun(@(x)isa(content(x,1),'Gridmod'),cconf);
            g = cconf{isg};
            gmod = g.content(1);
            
            self.K = cell(gmod.gridsize());
            self.H = cell(gmod.gridsize());
            self.gamma = cell(gmod.gridsize());
            self.mu = cell(gmod.gridsize());
            self.performance = cell(gmod.gridsize());
            self.solved = cell(gmod.gridsize());
            self.info = cell(gmod.gridsize());
            
            for k = 1:prod(gmod.gridsize())
                g.empty();
                g.add(gmod.grid_{k});
                [self.K{k},self.H{k},self.gamma{k},self.mu{k},self.performance{k},self.solved{k},self.info{k}] ...
                    = single_solve(self.cell_solver,config,specs,vars);
            end    
            self.K = Gridmod(self.K,gmod.params_);
            g.empty();
            g.add(gmod);
        end
    end
    
    methods (Static)
        function capabilities()
        % Dummy implementation for the abstract method of Solver. Has no
        % meaning for Gridsolver. 
            error('Gridmod is not a direct solver interface and, thus, has no capabilities on itself.');
        end
    end
end

