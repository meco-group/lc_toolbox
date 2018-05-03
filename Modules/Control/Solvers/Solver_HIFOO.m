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

classdef Solver_HIFOO < Solver
    % Solver_HIFOO defines a solver interface for \c hifoo.
    % \c hifoo solves the following problem:
    %
    % \image html gp_mixedFixedOrder.svg 
    %
    % @f[
    % \begin{aligned}
    % & \underset{K}{\text{minimize}} & & \max(\alpha_i\gamma_i, \ldots, \alpha_j\mu_j, \ldots) \\\
    % & \text{subject to} & & \left\lVert T\left(\frac{z_i}{w_i}\right)\right\rVert_\infty \leq \gamma_i, \quad i \text{ is an } \mathcal{H}_\infty \text{ objective}\\\
    % & & & \left\lVert T\left(\frac{z_j}{w_j}\right)\right\rVert_2 \leq \mu_i, \quad j \text{ is an } \mathcal{H}_2 \text{ objective} \\\
    % & & & \alpha_k \left\lVert T\left(\frac{z_k}{w_k}\right)\right\rVert_\infty \leq \beta_k, \quad k \text{ is an } \mathcal{H}_\infty \text{ constraint} \\\
    % & & & \alpha_l \left\lVert T\left(\frac{z_l}{w_l}\right)\right\rVert_2 \leq \beta_l, \quad l \text{ is an } \mathcal{H}_2 \text{ constraint} \\\
    % & & & T \text{ is stable} \\\
    % & & & K \text{ has order } p
    % \end{aligned}
    % @f]
    %
    % @f$P@f$ should be LTI, proper and stabilizable by @f$K@f$. 
    %
    % See \c hifoo.m for the syntax and for more details.
        
    methods
        function self = Solver_HIFOO(options)
        % Constructor for Solver_HIFOO objects.
        % 
        % Parameters:
        %  options: options that the user wants to pass to the solver
        %
        % Return values:
        %  self: the solver interface @type Solver_HIFOO 
        
            % Solver
            fprintf('Solving with HIFOO...\n\n');
            
            % Default options
            % none at this point

            if nargin > 0
                self = setoptions(self,options);
            end
        end
        
        function self = solve(self,config,specs,vars)
        % Parses the control problem such that it can be interpreted by \c
        % hifoo and calls the solver. 
        % 
        % Parameters:
        %  self: the solver interface @type Solver_HIFOO
        %  config: the control configuration @type SystemOfSystems
        %  specs: specifications of the control problem @type
        %  ControllerDesign
        %  vars: \c cell containing the optimization variables @type cell
        % 
        % Return values: 
        %  self: Solver object (hopefully) containing a solution to the
        %  control problem @type Solver_HIFOO
        
            ncont = length(specs.ctrl_in);
            nmeas = length(specs.ctrl_out);
            
            % HIFOO considers every channel as a separate generalized
            % plant, therefore split everything up
            P = cell(length(specs.performance),1);
            FUN = []; UPPERBND = []; OPTIONS = struct('weight',[]);
            for k = 1:length(specs.performance)
                thisspec = specs;
                thisspec.performance = thisspec.performance(k);
                thisplant = std(Solver.plant(config,thisspec,vars));
                thisplant.InputGroup.U1 = 1:length(thisspec.performance{:}.ch_p.in);
                thisplant.InputGroup.U2 = size(thisplant,2)-ncont+1:size(thisplant,2);
                thisplant.OutputGroup.Y1 = 1:length(thisspec.performance{:}.ch_p.out);
                thisplant.OutputGroup.Y2 = size(thisplant,1)-nmeas+1:size(thisplant,1);
                P{k} = thisplant;
                if isH2(thisspec.performance{:}); f = 't'; else; f = 'h'; end
                if k <= specs.nobj; u = Inf; else; u = upperbound(thisspec.performance{:}); end
                FUN = [FUN f];
                UPPERBND = [UPPERBND u];
                OPTIONS.weight = [OPTIONS.weight scale(specs.performance{k})];
            end
            
            % Controller order
            if specs.order == -1; GPfull = Solver.plant(config,specs,vars); ORDER = GPfull.content(1).nx; else; ORDER = specs.order; end
            
            % Solve
            tic
            [K, f, viol, loc] = hifoo(P, ORDER, FUN, UPPERBND, OPTIONS);
            self.info.time = toc;
            self.info.constrviol = viol;
            self.info.localoptimcert = loc;
            self.K = fromstd(K);
            
            warning('Since the objective is non-smooth, gamma and mu only consist of the worst-case norm that is present in the objective (can be H2 or H-infinity).');
            self.gamma = f;
            self.mu = f;
            self.solved = true;
            
            % Rescaling the performance only makes sense for the
            % constraints here
            specs = specs.rescale('constr+bound');
            
            % Throw warning
            if specs.nobj > 1
                warning('off','backtrace');
                warning('The objective for a multi-objective problem is formulated by HIFOO as obj = max(a*|.|, b*|.|, ...) and not as obj = a*|.| + b*|.| + ... (the way you formulated it). See the HIFOO documentation for more information.');
                warning('on','backtrace');
            end
            
        end
    end
    
    methods (Static)
        function cap = capabilities()
        % Returns the capabilities of \c HIFOO. 
        % 
        % Return values: 
        %  cap: capabilities of \c HIFOO @type struct
        
            cap.inout = 2;
            cap.norm = Inf;
            cap.constraints = true;
            cap.unstable = false;
            cap.improper = false;
            cap.parametric = false;
            cap.fixedorder = true;
        end
    end
end

