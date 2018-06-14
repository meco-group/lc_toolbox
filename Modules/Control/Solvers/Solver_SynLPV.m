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

classdef Solver_SynLPV < Solver
    % Solver_SynLPV defines a solver interface for \c LPV_unstructured_OF.
    % \c LPV_unstructured_OF solves one of the following problems:
    %
    % \image html gp_SynLPV.svg 
    %
    % @f[
    % \begin{aligned}
    % & \underset{K(\mathbf{p})}{\text{minimize}} & & \gamma \\\
    % & \text{subject to} & & \left\lVert T\left(\frac{z}{w}\right)\right\rVert_\infty \leq \gamma, \quad \forall \mathbf{p} \in \mathcal{P} \\\
    % & & & T \text{ is stable for a user-defined set of trajectories } \mathcal{P}
    % \end{aligned}
    % @f]
    %
    % @f[
    % \begin{aligned}
    % & \underset{K}{\text{minimize}} & & \mu \\\
    % & \text{subject to} & & \left\lVert T\left(\frac{z}{w}\right)\right\rVert_2 \leq \mu, \quad \forall \mathbf{p} \in \mathcal{P} \\\
    % & & & T \text{ is stable for a user-defined set of trajectories } \mathcal{P}
    % \end{aligned}
    % @f]
    %
    % @f$P(\mathbf{p})@f$ should be proper and stabilizable by @f$K(\mathbf{p})@f$ for the set of all allowed parameter trajectories @f$\mathcal{P}@f$. Its state-space matrices 
    % are allowed to be B-spline dependent on multiple scheduling parameters @f$\mathbf{p} = (p_1, p_2, \ldots)@f$. 
    % @f$\mathcal{P}@f$ is defined by the user and is characterized by bounds on the parameter @f$\mathbf{p}@f$ 
    % and its rate of variation @f$\mathbf{\dot p}@f$. The solution is found by transforming the problem into an
    % SDP containing sufficient LMI conditions.  
    % 
    % See \c 3rdparty/LCT_MixSynMIMO_LPV/LPV_unstructured_OF.m for the syntax and for more details.
    
    methods
        function self = Solver_SynLPV(options)
        % Constructor for Solver_SynLPV objects.
        % 
        % Parameters:
        %  options: options that the user wants to pass to the solver
        %
        % Return values:
        %  self: the solver interface @type Solver_SynLPV
        
            % Solver
            disp('Solving with SynLPV...');
            
            if nargin > 0
                self = setoptions(self,options);
            end
        end
        
        function self = solve(self,config,specs,vars)
        % Parses the control problem such that it can be interpreted by \c
        % LPV_unstructured_OF and calls the solver. 
        % 
        % Parameters:
        %  self: the solver interface @type Solver_SynLPV
        %  config: the control configuration @type SystemOfSystems
        %  specs: specifications of the control problem @type
        %  ControllerDesign
        %  vars: \c cell containing the optimization variables @type cell
        % 
        % Return values: 
        %  self: Solver object (hopefully) containing a solution to the
        %  control problem @type Solver_SynLPV
        
            % Compute plant state-space
            specs = rescale(specs,'all');
            P = Solver.plant(config,specs,vars);
            P = simplify(P);
            
            % solver input
            ncont = length(specs.ctrl_in);
            nmeas = length(specs.ctrl_out);
            if isH2(specs.performance{1})
                self.options.spec = 2;
            else
                self.options.spec = Inf;
            end
            [self.K,self.info] = LPV_unstructured_OF(P,nmeas,ncont,P.parameters(),self.options);
            
            % Save solver output
            self.gamma = self.info.objective;
            self.solved = true;
            self.performance{1} = specs.performance{1}*(1/self.gamma);
        end
%         
%         function bodemag(self)
%             error('Not implemented for LPV');
%         end
        
        function sigma(self)
            error('Not implemented for LPV');
        end
    end
    
    methods (Static)
        function cap = capabilities()
        % Returns the capabilities of \c SynLPV. 
        % 
        % Return values: 
        %  cap: capabilities of \c SynLPV @type struct
            cap.inout = 0;
            cap.norm = [2;Inf];
            cap.constraints = false;
            cap.unstable = false;
            cap.improper = false;
            cap.parametric = true;
            cap.fixedorder = false;
        end
    end
end