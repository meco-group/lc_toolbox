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
    %SOLVER_MIXEDFIXEDORDER Implementation of abstract solver for
    %mixedFixedOrder
    %   Solves both full and reduced order problems.
    
    methods
        function self = Solver_mixedFixedOrder(options)
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
            alpha(1,1:specs.nobj) = arrayfun(@(x) scale(x),specs.performance(1:specs.nobj));
            beta = zeros(1,length(specs.performance)); 
            beta(1,specs.nobj+1:end) = 1./arrayfun(@(x) scale(x),specs.performance(specs.nobj+1:end));
            
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
                info = compute_FO_mix(A,Bw,Bu,Cz,Cy,Dzw,Dzu,Dyw,Dyu,alpha,beta,ch,self.options)
                controllerfield = 'K';
                disp('SolveMixed full order controller computation')
            else
                FOspecs = specs;
                FOspecs.order = -1;
                FOsolver = Solver.select(config,FOspecs,vars,self.options);
                FOsolver = FOsolver.solve(config,FOspecs,vars);
                self.info.Kfo = fromstd(FOsolver.K);
                info = compute_RO_mix(A,Bw,Bu,Cz,Cy,Dzw,Dzu,Dyw,Dyu,...
                       balreal(std(self.info.Kfo)),specs.order,alpha,beta,ch,self.options,self.options)
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
                specs = specs.rescale('constr');
                self.performance = specs.performance;
                if specs.nobj > 0
                    obj = 1:specs.nobj;
                    gammaormu = self.gamma + self.mu;
                    self.performance(obj) = dealscale(self.performance(obj),(1./gammaormu(obj)));
                end
            else
                fprintf('Error while solving for the controller.\n');
                self.solved = false;
            end
        end
    end
    
    methods (Static)        
        function cap = capabilities()
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

