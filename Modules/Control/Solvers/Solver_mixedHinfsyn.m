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
    %SOLVER_MIXEDHINFSYN solve problem with mixedhinfsyn
    %   Abstract solver implementation for the mixedHinfsyn problem
        
    methods
        function self = Solver_mixedHinfsyn(options)
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
            % Compute plant state-space
            specs = specs.rescale('constr');
            P = Solver.plant(config,specs,vars);
            
            % Information to the solver
            alpha = zeros(length(specs.performance),1); % Setup which channels are objectives and which are constraints
            if(specs.nobj > 0) % Check if the problem is in fact an optimization problem
                alpha(1:specs.nobj,1) = arrayfun(@(x) scale(x),specs.performance(1:specs.nobj));
            end
            
            ncont = length(specs.ctrl_in);
            nmeas = length(specs.ctrl_out);
            P = std(P);
            if all(specs.performance.isoutput())
                Mz = arrayfun(@(x) size(x.W_out,1),specs.performance);
            else
                Mz = arrayfun(@(x) size(x.W_in,2),specs.performance);
                P = transpose(P);
            end
            Mz = transpose(Mz);
            
            % Compute the controller
            tic;
            [K,H,gamma] = mixedHinfsyn(balreal(P),nmeas,ncont,Mz,alpha,self.options);
            self.info.time = toc;
            K = balreal(K); % improve conditioning of sys_K
            
            % rescale performance weights
            self.performance = specs.performance;
            if specs.nobj > 0
                obj = 1:specs.nobj;
                self.performance(obj) = dealscale(self.performance(obj),1./gamma(obj,1));
            end
            
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