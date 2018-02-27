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

classdef Solver_Hinfsyn < Solver
    %SOLVER_HINFSYN Implementation of abstract solver class for hinfsyn
    %   Implementation of abstract solver class for hinfsyn
        
    methods
        function self = Solver_Hinfsyn(options)
            % Constructor for Solver_Hinfsyn objects
            %
            % Parameters: 
            % options : may have external options provided to the solver
            %
            % Return values:
            % self : options
            
            self.options.METHOD = 'lmi';
            
            if nargin > 0
                self = setoptions(self,options);
            end
        end
        
        function self = solve(self,config,specs,vars)
            % Compute plant state-space
            specs = rescale(specs,'all');
            P = Solver.plant(config,specs,vars);
            
            % solver input
            ncont = length(specs.ctrl_in);
            nmeas = length(specs.ctrl_out);

            opts = hinfsynopts(self);
            [K,H,gamma,self.info] = hinfsyn(balreal(std(P)),nmeas,ncont,opts{:});
            K = balreal(K); % improve conditioning of sys_K
            
            % rescale performance weights
            self.performance = specs.performance.*(1./gamma);
            
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
            cap.inout = 0;
            cap.norm = Inf;
            cap.constraints = false;
            cap.unstable = false;
            cap.improper = false;
            cap.parametric = false;
            cap.fixedorder = false;
        end
    end
    
    methods (Access=private)
        function list = hinfsynopts(self)
            list = {};
            keys = {'METHOD','GMAX','GMIN','TOLGAM','S0','DISPLAY'};
            for k = 1:length(keys)
                key = keys{k};
                if isfield(self.options,key) 
                    list = [list,{key self.options.(key)}];
                end
            end
        end
    end
end

