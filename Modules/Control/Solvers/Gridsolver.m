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

classdef Gridsolver < Solver
    %GRIDSOLVER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        cell_solver
    end
    
    methods
        function self = Gridsolver(cell_solver,options)
            self.cell_solver = cell_solver;
            self.options = options;
        end
        
        function self = solve(self,config,specs,vars)
            function [K,H,gamma,mu,performance,solved,info] = single_solve(solver,model,specs,vars)
                solution = solve(solver,model,specs,vars);
                K = solution.K;
                H = solution.H;
                gamma = solution.gamma;
                mu = solution.mu;
                performance = solution.performance;
                solved = solution.solved;
                info = solution.info;
            end
            
            cconf = content(remove_empty(unpack(config)));
            isg = cellfun(@(x)isa(content(x,1),'Gridmod'),cconf);
            g = cconf{isg};
            gmod = g.content(1);
            
            self.K = cell(gmod.gridsize());
            self.H = cell(gmod.gridsize());
            self.gamma = cell(gmod.gridsize());
            self.mu = cell(gmod.gridsize());
            self.performance = cell(gmod.gridsize());
            self.solved = cell(gmod.gridsize());
            self.info = cell(gmod.gridsize());
            
            for k = 1:prod(gmod.gridsize())
                g.empty();
                g.add(gmod.grid_{k});
                [self.K{k},self.H{k},self.gamma{k},self.mu{k},self.performance{k},self.solved{k},self.info{k}] ...
                    = single_solve(self.cell_solver,config,specs,vars);
            end    
            self.K = Gridmod(self.K,gmod.params_);
            g.empty();
            g.add(gmod);
        end
        
        function list = expand(self)
            function entry = makeentry(solver,k,h,h0,g,m,p,s,i,o)
                entry = solver;
                entry.K = k;
                
                c = SystemOfModels(h0.in,h0.out);
                c.add(h);
                entry.H = c;
                
                entry.gamma = g;
                entry.mu = m;
                entry.performance = p;
                entry.solved = s;
                entry.info = i;
                entry.options = o;
            end
            
            list = cellfun(@(k,h,g,m,p,s,i) makeentry(self.cell_solver,k,h,self.H,g,m,p,s,i,self.options), ...
                       self.K.grid_(:), self.H.content(1).grid_(:), self.gamma(:), self.mu(:), ...
                       self.performance(:), self.solved(:), self.info(:), 'un',0);
        end
    end
    
    methods (Static)
        function cap = capabilities()
            error('this function should not be called for a gridsolver');
        end
    end
end

