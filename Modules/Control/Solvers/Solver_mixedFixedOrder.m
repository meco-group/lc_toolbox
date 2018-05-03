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

classdef Solver_mixedFixedOrder < Solver
    % Solver_mixedFixedOrder defines a solver interface for \c compute_FO_mix and \c compute_RO_mix.
    % \c compute_FO_mix and \c compute_RO_mix together solve the following problem:
    %
    % \image html gp_mixedFixedOrder.svg 
    %
    % @f[
    % \begin{aligned}
    % & \underset{K}{\text{minimize}} & & \sum_i \alpha_i\gamma_i + \sum_j \alpha_j\mu_i \\\
    % & \text{subject to} & & \left\lVert T\left(\frac{z_i}{w_i}\right)\right\rVert_\infty \leq \gamma_i, \quad i \text{ is an } \mathcal{H}_\infty \text{ objective}\\\
    % & & & \left\lVert T\left(\frac{z_j}{w_j}\right)\right\rVert_2 \leq \mu_j, \quad j \text{ is an } \mathcal{H}_2 \text{ objective} \\\
    % & & & \left\lVert T\left(\frac{z_k}{w_k}\right)\right\rVert_\infty \leq 1, \quad k \text{ is an } \mathcal{H}_\infty \text{ constraint} \\\
    % & & & \left\lVert T\left(\frac{z_l}{w_l}\right)\right\rVert_2 \leq 1, \quad l \text{ is an } \mathcal{H}_2 \text{ constraint} \\\
    % & & & T \text{ is stable} \\\
    % & & & K \text{ has order } p
    % \end{aligned}
    % @f]
    %
    % @f$P@f$ should be LTI, proper and stabilizable by @f$K@f$. The
    % solution is found by transforming the problem into an SDP containing
    % sufficient LMI conditions. 
    %
    % See \c 3rdparty/LCT_MixedFixedOrder/compute_FO_mix.m and \c 3rdparty/LCT_MixedFixedOrder/compute_RO_mix.m for the syntax and for more details.
    
    methods
        function self = Solver_mixedFixedOrder(options)
        % Constructor for Solver_mixedFixedOrder objects.
        % 
        % Parameters:
        %  options: options that the user wants to pass to the solver
        %
        % Return values:
        %  self: the solver interface @type Solver_mixedFixedOrder
        
            % Solver
            disp('Solving with mixedFixedOrder')
            
            %Default options
            self.options.Dc = 0;
            self.options.Ts = 0;
            self.options.method = 'P';
            self.options.solver = 'mosek';
            self.options.verbose = 1;
            self.options.scaling = 1;
            
            if nargin > 0
                self = setoptions(self,options);
            end
        end
        
        function self = solve(self,config,specs,vars)
        % Parses the control problem such that it can be interpreted by \c
        % compute_FO_mix and \c compute_RO_mix and calls the solver(s). 
        % 
        % Parameters:
        %  self: the solver interface @type Solver_mixedFixedOrder
        %  config: the control configuration @type SystemOfSystems
        %  specs: specifications of the control problem @type
        %  ControllerDesign
        %  vars: \c cell containing the optimization variables @type cell
        % 
        % Return values: 
        %  self: Solver object (hopefully) containing a solution to the
        %  control problem @type Solver_mixedFixedOrder
        
            % Compute plant state-space
            % specs = specs.rescale('none');
            [P,wspecs] = Solver.plant(config,specs,vars);
            ch = Solver.channels(wspecs);
            
            P = balreal(std(P)); 
            if isdt(P)
                self.options.Ts = P.Ts;
                self.options.method = 'G';
            end
            
            % Provide input for solver
            alpha = zeros(1,length(specs.performance)); 
            alpha(1,1:specs.nobj) = cell2mat(cellfun(@(x) scale(x), specs.performance(1:specs.nobj), 'un', 0))';
            beta = zeros(1, length(specs.performance));
            beta(1,specs.nobj+1:end) = cell2mat(cellfun(@(x) upperbound(x)./scale(x), specs.performance(specs.nobj+1:end), 'un', 0))';
            
            % Split up the system in parts
            ncont = length(specs.ctrl_in);
            nmeas = length(specs.ctrl_out);
            nw = 1:(size(P,2)-ncont);
            nu = (nw(end)+1):size(P,2);
            nz = 1:(size(P,1)-nmeas);
            ny = (nz(end)+1):size(P,1);
            
            A = P.A; Bw = P.B(:,nw); Bu = P.B(:,nu);
            Cz = P.C(nz,:); Dzw = P.D(nz,nw); Dzu = P.D(nz,nu);
            Cy = P.C(ny,:); Dyw = P.D(ny,nw); Dyu = P.D(ny,nu);
            
            if specs.order == -1 % Solve full order problem
                info = compute_FO_mix(A,Bw,Bu,Cz,Cy,Dzw,Dzu,Dyw,Dyu,alpha,beta,ch,self.options);
                controllerfield = 'K';
                disp('SolveMixed full order controller computation')
            else
                FOspecs = specs;
                FOspecs.order = -1;
                FOsolver = Solver.select(config,FOspecs,vars,self.options);
                FOsolver = FOsolver.solve(config,FOspecs,vars);
                self.info.Kfo = fromstd(FOsolver.K);
                info = compute_RO_mix(A,Bw,Bu,Cz,Cy,Dzw,Dzu,Dyw,Dyu,...
                       balreal(std(self.info.Kfo)),specs.order,alpha,beta,ch,self.options,self.options);
                controllerfield = 'Kq';
                disp('SolveMixed reduced order controller computation')
            end    
            
            

            if info.feas                
                self.K = fromstd(info.(controllerfield));
                self.gamma = zeros(1,length(ch.In));
                self.mu = zeros(1,length(ch.In));
                if isfield(info,'gam')
                    self.gamma(1:length(info.gam)) = info.gam;
                end
                if isfield(info,'mu')
                    self.mu(1:length(info.mu)) = info.mu;
                end
                self.solved = true;
                self.info = mergestruct(self.info,info);
                
                % rescale performance
                specs = specs.rescale('constr+bound');
                self.performance = specs.performance;
                if specs.nobj > 0
                    obj = 1:specs.nobj;
                    gammaormu = self.gamma + self.mu;
                    self.performance(obj) = Norm.dealscale(self.performance(obj),num2cell(1./gammaormu(obj)'));
                end
            else
                self.solved = false;
                error('mixedFixedOrder was not able to find a feasible solution.');
            end
        end
    end
    
    methods (Static)        
        function cap = capabilities()
        % Returns the capabilities of \c mixedFixedOrder. 
        % 
        % Return values: 
        %  cap: capabilities of \c mixedFixedOrder @type struct
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

