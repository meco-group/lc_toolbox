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

classdef Solver_mixedHinfsyn < Solver
    % Solver_mixedHinfsyn defines a solver interface for \c mixedHinfsyn.
    % \c mixedHinfsyn solves the following problem:
    %
    % \image html gp_mixedHinfsyn.svg 
    %
    % @f[
    % \begin{aligned}
    % & \underset{K}{\text{minimize}} & & \sum_i \alpha_i\gamma_i^2 \\\
    % & \text{subject to} & & \left\lVert T\left(\frac{z_i}{w}\right)\right\rVert_\infty \leq \gamma_i, \quad i \text{ is an } \mathcal{H}_\infty \text{ objective} \\\
    % & & & \left\lVert T\left(\frac{z_j}{w}\right)\right\rVert_\infty \leq 1, \quad j \text{ is an } \mathcal{H}_\infty \text{ constraint} \\\
    % & & & T \text{ is stable}
    % \end{aligned}
    % @f]
    %
    % @f$P@f$ should be LTI, proper and stabilizable by @f$K@f$. The
    % solution is found by transforming the problem into an SDP containing
    % sufficient LMI conditions. 
    % 
    % See \c 3rdparty/LCT_MixedHinfSyn/mixedHinfsyn.m for the syntax and for more details.
        
    methods
        function self = Solver_mixedHinfsyn(options)
        % Constructor for Solver_mixedHinfsyn objects.
        % 
        % Parameters:
        %  options: options that the user wants to pass to the solver
        %
        % Return values:
        %  self: the solver interface @type Solver_mixedHinfsyn
        
            % Set default options
            self.options.gamma.solver = 'lmilab';
            self.options.gamma.solution = 2;
            self.options.controller.solution = 1;
            self.options.controller.LMItype = 1;
            self.options.controller.LMIsol = 2;
            
            if nargin > 0
                self = setoptions(self,options);
            end
        end
                
        function self = solve(self,config,specs,vars)
        % Parses the control problem such that it can be interpreted by \c
        % mixedHinfsyn and calls the solver. 
        % 
        % Parameters:
        %  self: the solver interface @type Solver_mixedHinfsyn
        %  config: the control configuration @type SystemOfSystems
        %  specs: specifications of the control problem @type
        %  ControllerDesign
        %  vars: \c cell containing the optimization variables @type cell
        % 
        % Return values: 
        %  self: Solver object (hopefully) containing a solution to the
        %  control problem @type Solver_mixedHinfsyn
            
            % Compute plant state-space
            specs = specs.rescale('constr+bound');
            P = Solver.plant(config,specs,vars);
            
            % Information to the solver
            alpha = zeros(length(specs.performance),1); % Setup which channels are objectives and which are constraints
            if(specs.nobj > 0) % Check if the problem is in fact an optimization problem
                alpha(1:specs.nobj,1) = cellfun(@(x) scale(x),specs.performance(1:specs.nobj));
            end
            
            ncont = length(specs.ctrl_in);
            nmeas = length(specs.ctrl_out);
            P = std(P);
            if all(cellfun(@isoutput, specs.performance))
                Mz = cellfun(@(x) size(x.W_out,1),specs.performance);
            else
                Mz = cellfun(@(x) size(x.W_in,2),specs.performance);
                P = transpose(P);
            end
            
            % Compute the controller
            tic;
            [K,H,gamma] = mixedHinfsyn(balreal(P),nmeas,ncont,Mz,alpha,self.options);
            self.info.time = toc;
            K = balreal(K); % improve conditioning of sys_K
            
            % Rescale performance weights
            self.performance = specs.performance;
            if specs.nobj > 0
                obj = 1:specs.nobj;
                self.performance(obj) = Norm.dealscale(self.performance(obj),num2cell(1./gamma(obj,1)));
            end
            
            % Save solver output
            self.K = fromstd(K);
            self.gamma = gamma;
            self.mu = zeros(size(gamma));
            self.solved = true;
            self.info.H = H;
            
            % Throw warning
            if specs.nobj > 0
                warning('off','backtrace');
                warning('mixedHinfsyn considers the squared norms as objective, i.e. obj = a*|.|² +  b*|.|² + ..., which is slightly different than the way you formulated it (without squares!). See the mixedHinfsyn documentation for more information.');
                warning('on','backtrace');
            end
        end
    end
    
    methods (Static)
        function cap = capabilities()
        % Returns the capabilities of \c mixedHinfsyn. 
        % 
        % Return values: 
        %  cap: capabilities of \c mixedHinfsyn @type struct
            cap.inout = 1;
            cap.norm = Inf;
            cap.constraints = true;
            cap.unstable = false;
            cap.improper = false;
            cap.parametric = false;
            cap.fixedorder = false;
        end
    end
end