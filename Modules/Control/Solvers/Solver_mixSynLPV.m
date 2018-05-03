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

classdef Solver_mixSynLPV < Solver
    % Solver_mixSynLPV defines a solver interface for \c LPV_unstructured_mix.
    % \c LPV_unstructured_mix solves the following problem:
    %
    % \image html gp_mixSynLPV.svg 
    %
    % @f[
    % \begin{aligned}
    % & \underset{K(\mathbf{p})}{\text{minimize}} & & \sum_i \alpha_i\gamma_i + \sum_j \alpha_j\mu_i \\\
    % & \text{subject to} & & \left\lVert T\left(\frac{z_i}{w_i}\right)\right\rVert_\infty \leq \gamma_i, \quad \forall \mathbf{p} \in \mathcal{P}, \quad i \text{ is an } \mathcal{H}_\infty \text{ objective}\\\
    % & & & \left\lVert T\left(\frac{z_j}{w_j}\right)\right\rVert_2 \leq \mu_i, \quad \forall \mathbf{p} \in \mathcal{P}, \quad j \text{ is an } \mathcal{H}_2 \text{ objective} \\\
    % & & & \left\lVert T\left(\frac{z_k}{w_k}\right)\right\rVert_\infty \leq 1, \quad \forall \mathbf{p} \in \mathcal{P}, \quad k \text{ is an } \mathcal{H}_\infty \text{ constraint} \\\
    % & & & \left\lVert T\left(\frac{z_l}{w_l}\right)\right\rVert_2 \leq 1, \quad \forall \mathbf{p} \in \mathcal{P}, \quad l \text{ is an } \mathcal{H}_2 \text{ constraint} \\\
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
    % See \c 3rdparty/LCT_MixSynMIMO_LPV/LPV_unstructured_mix.m for the syntax and for more details.
    
    methods
        function self = Solver_mixSynLPV(options)
        % Constructor for Solver_mixSynLPV objects.
        % 
        % Parameters:
        %  options: options that the user wants to pass to the solver
        %
        % Return values:
        %  self: the solver interface @type Solver_mixSynLPV
        
            % Solver
            disp('Solving with mixed SynLPV');
            self.options.objective = 'wc';
            
            if nargin > 0
                self = setoptions(self,options);
            end
        end
        
        function self = solve(self,config,specs,vars)
        % Parses the control problem such that it can be interpreted by \c
        % LPV_unstructured_mix and calls the solver. 
        % 
        % Parameters:
        %  self: the solver interface @type Solver_mixSynLPV
        %  config: the control configuration @type SystemOfSystems
        %  specs: specifications of the control problem @type
        %  ControllerDesign
        %  vars: \c cell containing the optimization variables @type cell
        % 
        % Return values: 
        %  self: Solver object (hopefully) containing a solution to the
        %  control problem @type Solver_mixSynLPV
        
            % Compute plant state-space
            specs = rescale(specs,'constr+bound');
            P = Solver.plant(config,specs,vars);
          
            % Weights for the norms in the objective
            alpha = zeros(1,length(specs.performance)); 
            alpha(1,1:specs.nobj) = cell2mat(cellfun(@(x) scale(x), specs.performance(1:specs.nobj), 'un', 0))';
            
            % Build the channel structure
            n_p = length(specs.performance);
            in_find = P.in;
            out_find = P.out; 
            in_index = 0;
            out_index = 0;
            channels(n_p) = struct('performance',[],'In',[],'Out',[]);
            for j = 1:n_p
                channels(j).performance = normtype(specs.performance{j});
                if j > specs.nobj
                    thisnorm = specs.performance{j}.norm;
                else
                    thisnorm = specs.performance{j};
                end
                if  isoutput(thisnorm)
                    [~,in_] = ismember(thisnorm.ch_p.in,in_find);
                    Nout = size(thisnorm.W_out,1);
                    out_ = out_index + (1:Nout);
                    out_index = out_index + Nout;
                else
                    [~,out_] = ismember(thisnorm.ch_p.out,out_find);
                    Nin = size(thisnorm.W_out,2);
                    in_ = in_index + (1:Nin);
                    in_index = in_index + Nin;
                end
                channels(j).In = in_;
                channels(j).Out = out_;
            end

            % Solver input
            ncont = length(specs.ctrl_in);
            nmeas = length(specs.ctrl_out);
            P = simplify(P.content(1));
            [self.K,self.info] = LPV_unstructured_mix(P,nmeas,ncont,alpha,channels,P.parameters(),self.options);
            self.solved = true;
            
            % Save solver output
            self.performance = specs.performance;
            self.gamma = zeros(1,length(self.performance));
            self.mu = zeros(1,length(self.performance));
            
            for k = 1:length(channels)
                if k <= specs.nobj
                    val = sqrt(self.info.objective);
                else
                    val = 1;
                end
                
                if channels(k).performance == Inf
                    self.gamma(k) = val;
                else
                    self.mu(k) = val;
                end
            end
            
            if strcmp(self.options.objective,'wc')
                gamormu = self.gamma + self.mu;
                if specs.nobj > 0
                    obj = 1:specs.nobj;
                    self.performance(obj) = cellfun(@(x,y) {x*(1/y)}, self.performance(obj), num2cell(gamormu(1:specs.nobj)'));
                end
            end
        end
    end   
        methods (Static)
        function cap = capabilities()
        % Returns the capabilities of \c mixSynLPV. 
        % 
        % Return values: 
        %  cap: capabilities of \c mixSynLPV @type struct
            cap.inout = 2;
            cap.norm = [2;Inf];
            cap.constraints = true;
            cap.unstable = false;
            cap.improper = false;
            cap.parametric = true;
            cap.fixedorder = false;
        end
    end
end