% This file is part of LCToolbox.
% (c) Copyright 2024 - MECO Research Team, KU Leuven. 
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

classdef Solver_gssshinfcd < Solver
    
    properties
        solverobj; 
    end
    
    methods
        function self = Solver_gssshinfcd(options)
        
            % Solver
            disp('Solving with gssshinfcd'); 
            
            self.solverobj = gssshinfcd.gssshinfcd(); 
            self.options = self.solverobj.options(); 
            
            if nargin > 0
                self.options = mergestruct(options,self.options); 
            end
            
        end
        
        function self = solve(self,config,specs,vars)
        
            % Compute plant state-space
            specs = rescale(specs,'constr+bound');
            [P,~,~,ch] = Solver.plant(config,specs,vars);
            P = simplify(P); 
            assert(isa(P,'LPVDSSmod'), 'Could not transform your generalized plant model to a state-space representation with a B-spline dependency on the parameters.'); 
            ncont = length(specs.ctrl_in);
            nmeas = length(specs.ctrl_out);
          
            % Information to the solver
            P = gssshinfcd.gsss.fromlpvdss(P);
            CH = struct('in',[],'out',[],'weight',[]);
            for i=1:length(ch.Hinf)
                [CH(i).in,~] = find(ch.In{i});
                [~,CH(i).out] = find(ch.Out{i}); 
                if i<=specs.nobj
                    CH(i).weight = scale(specs.performance{i}); 
                else
                    CH(i).weight = 0;
                end
            end        

            % Configure the solver
            if isempty(specs.region)
                self.solverobj = self.solverobj.setproblem(P,nmeas,ncont,CH);
            else
                for i=1:length(specs.region)
                    reg(i).M = specs.region{i}.M;
                    reg(i).L = specs.region{i}.L; 
                end
                self.solverobj = self.solverobj.setproblem(P,nmeas,ncont,CH,reg);
            end
            self.solverobj = self.solverobj.setoptions(self.options);
            
            % Compute the controller
            tic;
            K = self.solverobj.solve();
            self.info.time = toc;
            
%            % Rescale performance weights % -> not possible with LPV... maybe in LFT form? -> to be investigated later
%             self.performance = specs.performance;
%             if specs.nobj > 0
%                 objectives = 1:specs.nobj;
%                 self.performance(objectives) = Norm.dealscale(self.performance(objectives),num2cell(1./sqgamma(objectives,1)));
%             end
            
            % Save solver output
            self.K = lft(K(1).K.tolpvdss,K(1).fb); 
            self.gamma = K(1).sqgamma;
            self.mu = [];
            self.solved = true;
            
        end
    end   
        methods (Static)
        function cap = capabilities()
        % Returns the capabilities of \c mixSynLPV. 
        % 
        % Return values: 
        %  cap: capabilities of \c mixSynLPV @type struct
            cap.inout = 2;
            cap.norm = [2;Inf];
            cap.constraints = true;
            cap.unstable = false;
            cap.improper = false;
            cap.parametric = true;
            cap.fixedorder = false;
            cap.polereg = false;
        end
    end
end