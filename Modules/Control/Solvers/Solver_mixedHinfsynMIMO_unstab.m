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

classdef Solver_mixedHinfsynMIMO_unstab < Solver
    %SOLVER_MIXEDHINFSYNMIMO_unstab Implementation of abstract solver class for
    %mixedHinfsynMIMO_unstab
    %   Detailed explanation goes here
        
    methods
        function self = Solver_mixedHinfsynMIMO_unstab(options)
            % Solver
            disp('Solving with mixedHinfsynMIMO_unstab')
            
            % Default options
            self.options.gammasolver = 'lmilab';
            self.options.controllersolver = 'lmilab';

            if nargin > 0
                self = setoptions(self,options);
            end
        end
        
        function self = solve(self,config,specs,vars)
            % Get output filters
            specs = rescale(specs,'constr+bound');
            [P,Wu,ch] = Solver_mixedHinfsynMIMO_unstab.plant(config,specs,vars);
    
            % Setup which channels are objectives and which are constraints
            alpha = zeros(1,length(specs.performance)); 
            alpha(1,1:specs.nobj) = cell2mat(cellfun(@(x) scale(x), specs.performance(1:specs.nobj), 'un', 0))';
            
            % Get number of controls and measurements
            ncont = length(specs.ctrl_in);
            nmeas = length(specs.ctrl_out);

            tic;
            if isnumeric(Wu); Wu = fromstd(Wu); end
            [K,gamma] = mixedHinfsyn_MIMO(balreal(std(P)),std(Wu),nmeas,ncont,alpha,ch,self.options);
            self.info.time = toc;
            self.K = fromstd(K);
            self.gamma = transpose(gamma(:,1)); 
            self.mu = zeros(size(gamma)); 
            self.solved = true;
            
            % rescale performance weights
            self.performance = specs.performance;
            if specs.nobj > 0
                obj = 1:specs.nobj;
                self.performance(obj) = Norm.dealscale(self.performance(obj),num2cell(1./self.gamma(obj)));
            end
        end
    end
    
    methods (Static)
        function [P,Wu,ch] = plant(config,specs,vars)
            % Separates stable and unstable dynamics in the weight
            
            Wu = [];
            stabspecs = specs;
            % set up unstable weights
            for k = 1:length(specs.performance)
                [GS,GNS] = stabsep(specs.performance{k}.W_out);
                if GNS.nx == 0
                    Wu = blkdiag(Wu,eye(size(specs.performance{k}.W_out,1)));
                elseif GS.nx == 0
                    Wu = blkdiag(Wu,specs.performance{k}.W_out);
                    stabspecs.performance{k} = setWout(stabspecs.performance{k},eye(size(specs.performance{k}.W_out,1)));
                else
                    error('Are you kidding me? I detected both stable and unstable poles in one weight.\n In case you want an integrator weight, make all poles unstable. Otherwise, make everything stable.')
                end
                assert(isstable(specs.performance{k}.W_in),'mixedHinfsynMIMO_unstab does not support unstable input weighting filters. Try mixedHinfsynMIMO instead.');
            end
                
            [P,wspecs] = Solver.plant(config,stabspecs,vars);
            ch = Solver.channels(wspecs);            
        end
        
        function cap = capabilities()
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