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

classdef Solver_systune < Solver
    % Solver_systune defines a solver interface for \c systune.
    % \c systune solves the following problem:
    %
    % \image html gp_mixedFixedOrder.svg 
    %
    % @f[
    % \begin{aligned}
    % & \underset{K}{\text{minimize}} & & \max(\gamma_i, \ldots,\mu_j, \ldots) \\\
    % & \text{subject to} & & \left\lVert T\left(\frac{z_i}{w_i}\right)\right\rVert_\infty \leq \gamma_i, \quad i \text{ is an } \mathcal{H}_\infty \text{ objective}\\\
    % & & & \left\lVert T\left(\frac{z_j}{w_j}\right)\right\rVert_2 \leq \mu_j, \quad j \text{ is an } \mathcal{H}_2 \text{ objective} \\\
    % & & & \left\lVert T\left(\frac{z_k}{w_k}\right)\right\rVert_\infty \leq \beta_k, \quad k \text{ is an } \mathcal{H}_\infty \text{ constraint} \\\
    % & & & \left\lVert T\left(\frac{z_l}{w_l}\right)\right\rVert_2 \leq \beta_l, \quad l \text{ is an } \mathcal{H}_2 \text{ constraint} \\\
    % & & & T \text{ is stable} \\\
    % & & & K \text{ has order } p
    % \end{aligned}
    % @f]
    %
    % @f$P@f$ should be LTI, proper and stabilizable by @f$K@f$. 
    %
    % See \c systune (MATLAB's Robust Control Toolbox) for the syntax and for more details.
        
    methods
        function self = Solver_systune(~)
        % Constructor for Solver_systune objects.
        % 
        % Return values:
        %  self: the solver interface @type Solver_systune
        
            % Do nothing for now - options are ignored.
        end
        
        function self = solve(self,config,specs,vars)
        % Parses the control problem such that it can be interpreted by \c
        % systune and calls the solver. 
        % 
        % Parameters:
        %  self: the solver interface @type Solver_systune
        %  config: the control configuration @type SystemOfSystems
        %  specs: specifications of the control problem @type
        %  ControllerDesign
        %  vars: \c cell containing the optimization variables @type cell
        % 
        % Return values: 
        %  self: Solver object (hopefully) containing a solution to the control problem @type Solver_systune
        
            % Compute plant state-space
            if ~(all(cellfun(@(x) scale(x) == 1, specs.performance(1:specs.nobj))))
                warning('Systune rescales the input weights by the provided scaling factors.');
            end
            specs = specs.rescale('obj');
            [P,wspecs] = Solver.plant(config,specs,vars);
            ch = Solver.channels(wspecs);   

            % Split up the system in parts
            ncont = length(specs.ctrl_in);
            nmeas = length(specs.ctrl_out);
            nw = 1:(size(P,2)-ncont);
            nu = (nw(end)+1):size(P,2);
            nz = 1:(size(P,1)-nmeas);
            ny = (nz(end)+1):size(P,1);

            % define systune tuning goals   
            P = std(P);
            P.inputname = arrayfun(@(x)['u' num2str(x)],1:size(P,2),'UniformOutput',false);
            P.outputname = arrayfun(@(x)['y' num2str(x)],1:size(P,1),'UniformOutput',false);
            objectives = []; constraints = [];
            for k = 1:length(ch.In)
                inputs = arrayfun(@(x)['u' num2str(x)],find(sum(ch.In{k},2)),'UniformOutput',false);
                outputs = arrayfun(@(x)['y' num2str(x)],find(sum(ch.Out{k},1)),'UniformOutput',false);
                
                if isH2(specs.performance{k}); tg = 'Variance'; else; tg = 'Gain'; end
                if k <= specs.nobj
                    sc = 1/scale(specs.performance{k});
                    objectives = [objectives TuningGoal.(tg)(inputs,outputs,sc)];
                else
                    sc = upperbound(specs.performance{k})/scale(specs.performance{k});
                    constraints = [constraints TuningGoal.(tg)(inputs,outputs,sc)];
                end
            end

            % make controller block
            if specs.order == -1
                K = ltiblock.ss('K',size(P.a,1),length(ny),length(nu));
            else
                K = ltiblock.ss('K',specs.order,length(ny),length(nu));
            end
            K.u = P.outputname(ny);
            K.y = P.inputname(nu);

            % Try solving for a controller
            P = connect(P,K,P.inputname(nw),P.outputname(nz));
            [CL,fSoft,gHard,info] = systune(P,objectives,constraints);
            gammaormu = [fSoft,gHard];
            
            % rescale performance
            self.performance = specs.performance;
            if specs.nobj > 0
                obj = 1:specs.nobj;
                self.performance(obj) = Norm.dealscale(self.performance(obj),num2cell(1./fSoft(obj)));
            end
            
            % save output
            self.K = fromstd(getBlockValue(CL,'K'));
            self.gamma = gammaormu;
            self.gamma(ch.H2) = 0;
            self.mu = gammaormu;
            self.mu(ch.Hinf) = 0;
            self.info = mergestruct(self.info,info);
            self.solved = true;
            
            % throw warning
            if specs.nobj > 1
                warning('off','backtrace');
                warning('The objective for a multi-objective problem is formulated by systune as obj = max(a*|.|, b*|.|, ...) and not as obj = a*|.| + b*|.| + ... (the way you formulated it). See the systune documentation for more information.');
                warning('on','backtrace');
            end
        end
    end
    
    methods (Static)
        function cap = capabilities()
        % Returns the capabilities of \c systune. 
        % 
        % Return values: 
        %  cap: capabilities of \c systune @type struct
            cap.inout = 2;
            cap.norm = [2;Inf];
            cap.constraints = true;
            cap.unstable = false;
            cap.improper = false;
            cap.parametric = false;
            cap.fixedorder = true;
        end
    end
end

