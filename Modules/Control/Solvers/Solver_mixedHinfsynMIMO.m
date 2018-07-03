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

classdef Solver_mixedHinfsynMIMO < Solver
    % Solver_mixedHinfsynMIMO defines a solver interface for \c mixedHinfsynMIMO.
    % \c mixedHinfsynMIMO solves the following problem:
    %
    % \image html gp_mixedHinfsynMIMO.svg 
    %
    % @f[
    % \begin{aligned}
    % & \underset{K}{\text{minimize}} & & \sum_i \alpha_i\gamma_i \\\
    % & \text{subject to} & & \left\lVert T\left(\frac{z_i}{w_i}\right)\right\rVert_\infty \leq \gamma_i, \quad i \text{ is an } \mathcal{H}_\infty \text{ objective}\\\\
    % & & & \left\lVert T\left(\frac{z_j}{w_j}\right)\right\rVert_\infty\leq \beta_j, \quad j \text{ is an } \mathcal{H}_\infty \text{ constraint}\ \\\
    % & & & T \text{ is input-output stable (not necessarily internally stable)}
    % \end{aligned}
    % @f]
    %
    % @f$P@f$ should be LTI, proper and stabilizable by @f$K@f$. This solver, 
    % however, allows @f$W_i@f$ and @f$W_o@f$ to contain unstable modes
    % (e.g. weights containing integrators) and reduces the controller
    % order without performance loss in case of a singular problem.  The
    % solution is found by transforming the problem into an SDP containing
    % sufficient LMI conditions. 
    %
    % See \c 3rdparty/LCT_MixedHinfSynMIMO/mixedHinfsynMIMO.m for the syntax and for more details.
        
    methods
        function self = Solver_mixedHinfsynMIMO(options)
        % Constructor for Solver_mixedHinfsynMIMO objects.
        % 
        % Parameters:
        %  options: options that the user wants to pass to the solver
        %
        % Return values:
        %  self: the solver interface @type Solver_mixedHinfsynMIMO 
        
            % Solver
            fprintf('Solving with mixedHinfsynMIMO...\n\n');
            
            % Default options
            self.options.gammasolver.solver = 'mosek';
            self.options.gammasolver.mosek.MSK_DPAR_ANA_SOL_INFEAS_TOL = 1e-8;
            self.options.controllersolver.solver = 'basiclmi';

            if nargin > 0
                self = setoptions(self,options);
            end
        end
        
        function self = solve(self,config,specs,vars)
        % Parses the control problem such that it can be interpreted by \c
        % mixedHinfsynMIMO and calls the solver. 
        % 
        % Parameters:
        %  self: the solver interface @type Solver_mixedHinfsynMIMO
        %  config: the control configuration @type SystemOfSystems
        %  specs: specifications of the control problem @type
        %  ControllerDesign
        %  vars: \c cell containing the optimization variables @type cell
        % 
        % Return values: 
        %  self: Solver object (hopefully) containing a solution to the
        %  control problem @type Solver_mixedHinfsynMIMO
        
            % Rescale constraints
            %specs = rescale(specs,'constr+bound');
            
            % Get output filters
            [P,Wo,Wi,ch] = Solver_mixedHinfsynMIMO.plant(config,specs,vars);
    
            % Setup which channels are objectives and which are constraints
            alpha = zeros(1,length(specs.performance)); 
            alpha(1,1:specs.nobj) = cell2mat(cellfun(@(x) scale(x), specs.performance(1:specs.nobj), 'un', 0))';
            beta = zeros(1, length(specs.performance));
            beta(1,specs.nobj+1:end) = cell2mat(cellfun(@(x) upperbound(x)./scale(x), specs.performance(specs.nobj+1:end), 'un', 0))';
            
            % Get number of controls and measurements
            ncont = length(specs.ctrl_in);
            nmeas = length(specs.ctrl_out);

            tic;
%             stdP = std(P);
%             if isnumeric(Wo); Wo = ss(Wo); Wo.Ts = stdP.Ts; else Wo = std(Wo); end
%             if isnumeric(Wi); Wi = ss(Wi); Wi.Ts = stdP.Ts; else Wi = std(Wi); end
            [K, gamma, ~] = mixedHinfsynMIMO(P,Wi,Wo,nmeas,ncont,alpha,beta,ch,self.options);
            self.info.time = toc;
            self.K = fromstd(K);
            self.gamma = gamma(:,1); 
            self.mu = zeros(size(gamma)); 
            self.solved = true;
            
            % rescale performance
            specs = rescale(specs,'constr+bound');
            self.performance = specs.performance;
            if specs.nobj > 0
                obj = 1:specs.nobj;
                self.performance(obj) = Norm.dealscale(self.performance(obj),num2cell(1./self.gamma(obj,1)));
            end
        end
    end
    
    methods (Static)
        function [P,Wo,Wi,ch] = plant(config,specs,vars)
        % Parses the generalized plant in the form that is required by
        % \c mixedHinfsynMIMO. 
        % 
        % Parameters:
        %  config: the control configuration @type SystemOfSystems
        %  specs: specifications of the control problem @type
        %  ControllerDesign
        %  vars: \c cell containing the optimization variables @type cell
        %
        % Return values:
        %  P: generalized plant @type numlti
        %  Wo: unstable output weighting filter @type numlti
        %  Wi: unstable input weighting filter @type numlti
        %  ch: structure defining the channels corresponding to the
        %  specifications @type struct
        
%             Wo = [];
%             Wi = [];
%             stabspecs = specs;
%             
%             % Set up unstable weights
%             for k = 1:length(specs.performance)
%                 
%                 [GSo,GNSo] = stabsep(fromstd(specs.performance{k}.W_out));
%                 [GSi,GNSi] = stabsep(fromstd(specs.performance{k}.W_in));
%                 if GNSo.nx == 0 % stable weight
%                     Wo = blkdiag(Wo,eye(size(specs.performance{k}.W_out,1)));
%                 elseif GSo.nx == 0 % unstable weight
%                     Wo = blkdiag(Wo,specs.performance{k}.W_out);
%                     stabspecs.performance{k} = setWout(stabspecs.performance{k},eye(size(specs.performance{k}.W_out,1)));
%                 else
%                     error('One of your output weights contains both stable and unstable poles, which I cannot separate. Make all poles unstable in case you want a weight with integrators and make all poles stable otherwise.')
%                 end
%                 if GNSi.nx == 0
%                     Wi = blkdiag(Wi,eye(size(specs.performance{k}.W_in,1)));
%                 elseif GSi.nx == 0
%                     Wi = blkdiag(Wi,specs.performance{k}.W_in);
%                     stabspecs.performance{k} = setWin(stabspecs.performance{k},eye(size(specs.performance{k}.W_in,1)));
%                 else
%                     error('One of your input weights contains both stable and unstable poles, which I cannot separate. Make all poles unstable in case you want a weight with integrators and make all poles stable otherwise.')
%                 end
%             end
%                 
            [P,Wo,Wi,ch] = Solver.plant(config,specs,vars);
            P = std(P);
%             ch = Solver.channels(wspecs,P);   
            for i = 1:length(ch.In)
                [r,~] = find(ch.In{i});
                [~,c] = find(ch.Out{i});
                channels(i).In = r';
                channels(i).Out = c';
            end
            ch = channels;
        end
        
        function cap = capabilities()
        % Returns the capabilities of \c mixedHinfsynMIMO. 
        % 
        % Return values: 
        %  cap: capabilities of \c mixedHinfsynMIMO @type struct
        
            cap.inout = 2;
            cap.norm = Inf;
            cap.constraints = true;
            cap.unstable = true;
            cap.improper = false;
            cap.parametric = false;
            cap.fixedorder = false;
        end
    end
end

