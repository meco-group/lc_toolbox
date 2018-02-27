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
    %SOLVER_SYNLPV Summary of this class goes here
    %   Detailed explanation goes here
    
    methods
        function self = Solver_SynLPV(options)
            % Solver
            disp('Solving with SynLPV');
            
            if nargin > 0
                self = setoptions(self,options);
            end
        end
        
        function self = solve(self,config,specs,vars)
            % Compute plant state-space
            specs = rescale(specs,'all');
            P = Solver.plant(config,specs,vars);
            P = simplify(P.content(1));
            
            % solver input
            ncont = length(specs.ctrl_in);
            nmeas = length(specs.ctrl_out);
            if isH2(specs.performance(1))
                self.options.spec = 2;
            else
                self.options.spec = Inf;
            end
            [self.K,self.info] = LPV_unstructured_OF(P,nmeas,ncont,P.parameters(),self.options);
            
            % Save solver output
            self.gamma = self.info.objective;
            self.solved = true;
            self.performance = specs.performance*(1/self.gamma);
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