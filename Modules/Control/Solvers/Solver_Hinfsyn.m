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

classdef Solver_Hinfsyn < Solver
    % Solver_Hinfsyn defines a solver interface for \c hinfsyn.
    % \c hinfsyn solves the following problem:
    %
    % \image html gp_Hinfsyn.svg 
    %
    % @f[
    % \begin{aligned}
    % & \underset{K}{\text{minimize}} & & \gamma \\\
    % & \text{subject to} & & \left\lVert T\left(\frac{z}{w}\right)\right\rVert_\infty \leq \gamma\\\
    % & & & T \text{ is stable}
    % \end{aligned}
    % @f]
    %
    % @f$P@f$ should be LTI, proper and stabilizable by @f$K@f$. 
    % 
    % See \c hinfsyn (MATLAB's Robust Control Toolbox) for the syntax and for more details.
        
    methods
        function self = Solver_Hinfsyn(options)
        % Constructor for Solver_Hinfsyn objects.
        % 
        % Parameters:
        %  options: options that the user wants to pass to the solver
        %
        % Return values:
        %  self: the solver interface @type Solver_Hinfsyn
            
            self.options.METHOD = 'lmi'; % default
            
            if nargin > 0
                self = setoptions(self,options);
            end
        end
        
        function self = solve(self,config,specs,vars)
        % Parses the control problem such that it can be interpreted by \c
        % hinfsyn and calls the solver. 
        % 
        % Parameters:
        %  self: the solver interface @type Solver_Hinfsyn
        %  config: the control configuration @type SystemOfSystems
        %  specs: specifications of the control problem @type
        %  ControllerDesign
        %  vars: \c cell containing the optimization variables @type cell
        % 
        % Return values: 
        %  self: Solver object (hopefully) containing a solution to the
        %  control problem @type Solver_Hinfsyn
        
            % Compute plant state-space
            specs = rescale(specs,'all');
            P = Solver.plant(config,specs,vars);
            
            % solver input
            ncont = length(specs.ctrl_in);
            nmeas = length(specs.ctrl_out);

            opts = hinfsynopts(self);
            [K,H,gamma,self.info] = hinfsyn(balreal(std(P)),nmeas,ncont,opts{:});
            K = balreal(K); % improve conditioning of sys_K
            
            % rescale performance weights
            self.performance = cellfun(@(x) {x*(1/gamma)}, specs.performance);
            
            % Save solver output
            self.K = fromstd(K);
            self.gamma = gamma;
            self.mu = zeros(size(gamma));
            self.solved = true;
            self.info.H = H;
        end
    end
    
    methods (Static)
        function cap = capabilities()
            cap.inout = 0;
            cap.norm = Inf;
            cap.constraints = false;
            cap.unstable = false;
            cap.improper = false;
            cap.parametric = false;
            cap.fixedorder = false;
        end
    end
    
    methods (Access=private)
        function list = hinfsynopts(self)
            list = {};
            keys = {'METHOD','GMAX','GMIN','TOLGAM','S0','DISPLAY'};
            for k = 1:length(keys)
                key = keys{k};
                if isfield(self.options,key) 
                    list = [list,{key self.options.(key)}];
                end
            end
        end
    end
end

