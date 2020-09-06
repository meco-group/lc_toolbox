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

classdef Solver_hinfcd < Solver
    
    properties
        solverobj; 
    end
        
    methods
        function self = Solver_hinfcd(options)
            disp('solving with hinfcd');
            
            self.solverobj = hinfcd.hinfcd(); 
            self.options = self.solverobj.options(); 
            
            if nargin > 0
                self.options = mergestruct(options,self.options); 
                self = setoptions(self,options);
            end
        end
                
        function self = solve(self,config,specs,vars)
            
            % Compute plant state-space
            specs = specs.rescale('constr+bound');
            opts.separate_filters = true; 
            [P,Wout,Win,ch] = Solver.plant(config,specs,vars,opts);
            
            % Information to the solver
            P = std(P);
            Wout = std(Wout);
            Win = std(Win);
            
            CH = struct('in',[],'out',[],'weight',[]);
            for i=1:length(ch.Hinf)
                CH(i).in = find(ch.In{i});
                CH(i).out = find(ch.Out{i}); 
                if i<=specs.nobj
                    CH(i).weight = scale(specs.performance{i}); 
                else
                    CH(i).weight = 0;
                end
            end
            
            % Configure the solver
            self.solverobj = self.solverobj.setproblem(P,Win,Wout,CH);
            
            % Compute the controller
            tic;
            [K,sqgamma] = self.solverobj.solve(self.options);
            self.info.time = toc;
            
            % Rescale performance weights
            self.performance = specs.performance;
            if specs.nobj > 0
                objectives = 1:specs.nobj;
                self.performance(objectives) = Norm.dealscale(self.performance(objectives),num2cell(1./sqgamma(objectives,1)));
            end
            
            % Save solver output
            self.K = fromstd(K);
            self.gamma = sqrt(sqgamma);
            self.mu = zeros(size(sqgamma));
            self.solved = true;
            
        end
    end
    
    methods (Static)
        function cap = capabilities()
            cap.inout = 2;
            cap.norm = Inf;
            cap.constraints = true;
            cap.unstable = true;
            cap.improper = true;
            cap.parametric = false;
            cap.fixedorder = false;
        end
    end
end
