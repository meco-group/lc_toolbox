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
    %SOLVER_MIXEDHINFSYNMIMO solve problem with mixedHinfsynMIMO
    %   Abstract solver implementation for the mixedHinfsynMIMO problem
        
    methods
        function self = Solver_mixedHinfsynMIMO(options)
            % Set default options
            self.options.gammasolver = 'mosek';
            self.options.controllersolver = 'basiclmi';
            
            if nargin > 0
                self = setoptions(self,options);
            end
        end
                
        function self = solve(self,config,specs,vars)  
            % Compute plant state-space
            specs = specs.rescale('constr');
            [P,ch] = self.plant(config,specs,vars);
            
            % Information to the solver
            alpha = zeros(length(specs.performance),1); % Setup which channels are objectives 
            if(specs.nobj > 0) % Check if the problem is in fact an optimization problem
                alpha(1:specs.nobj,1) = arrayfun(@(x) scale(x),specs.performance(1:specs.nobj));
            end
            
            beta = zeros(length(specs.performance),1); % Setup which channels are objectives and which are constraints
            if(length(specs.performance) - specs.nobj > 0) % Check if the problem has in fact constraints
                beta(specs.nobj+1:end,1) = arrayfun(@(x) scale(x),specs.performance(specs.nobj+1:end));
            end
            
            ncont = length(specs.ctrl_in);
            nmeas = length(specs.ctrl_out);
            P = balreal(std(P));
            
            % Compute the controller
            tic;
            [K,gamma] = mixedHinfsynMIMO(P,nmeas,ncont,alpha,beta,ch,self.options);
            self.info.time = toc;
            K = balreal(K); % improve conditioning of sys_K
            
            % rescale performance weights
            self.performance = specs.performance;
            if specs.nobj > 0
                obj = 1:specs.nobj;
                self.performance(obj) = dealscale(self.performance(obj),1./gamma(1,obj));
            end
            
            % Save solver output
            self.K = fromstd(K);
            self.gamma = gamma;
            self.mu = zeros(size(gamma));
            self.solved = true;
        end
    end
    
    methods (Static)
        
        function [P,ch] = plant(config,specs,vars)
            % translate into alternative ('new') channel representation
            [P,wspecs] = Solver.plant(config,specs,vars);
            ch = Solver.channels(wspecs);   
            for i = 1:length(ch.In)
                [r,~] = find(ch.In{i});
                [~,c] = find(ch.Out{i});
                channels(i).In = r';
                channels(i).Out = c';
            end
            ch = channels;
        end
        
        function cap = capabilities()
            cap.inout = 2;
            cap.norm = Inf;
            cap.constraints = true;
            cap.unstable = false;
            cap.improper = false;
            cap.parametric = false;
            cap.fixedorder = false;
        end
    end
end