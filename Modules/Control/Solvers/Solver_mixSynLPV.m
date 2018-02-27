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

classdef Solver_mixSynLPV < Solver
    %SOLVER_SYNLPV Summary of this class goes here
    %   Detailed explanation goes here
    
    methods
        function self = Solver_mixSynLPV(options)
            % Solver
            disp('Solving with mixed SynLPV');
            self.options.objective = 'wc';
            
            if nargin > 0
                self = setoptions(self,options);
            end
        end
        
        function self = solve(self,config,specs,vars)
            % Compute plant state-space
            specs = rescale(specs,'all');
            [P,wspecs] = Solver.plant(config,specs,vars);
            ch = Solver.channels(wspecs);
            % Information to solver
            alpha = zeros(1,length(specs.performance)); 
            alpha(1,1:specs.nobj) = 1;
            n_p = length(specs.performance);
            in_find = P.in;
            out_find = P.out; 
            in_index = 0;
            out_index = 0;
            channels(n_p) = struct('performance',[],'In',[],'Out',[]);
                for j = 1:n_p
                    channels(j).performance = specs.performance(j).p;
                    
                    if all(specs.performance(j).isoutput())
                        [~,in_] = ismember(specs.performance(j).ch_p.in,in_find);
                        Nout = size(specs.performance(j).W_out,1);
                        out_ = out_index + (1:Nout);
                        out_index = out_index + Nout;
                    else
                        [~,out_] = ismember(specs.performance(j).ch_p.out,out_find);
                        Nin = size(specs.performance(j).W_out,2);
                        in_ = in_index + (1:Nin);
                        in_index = in_index + Nin;
                    end
                    channels(j).In = in_;
                    channels(j).Out = out_;
                end
                        %channels = ch;
            % solver input
            ncont = length(specs.ctrl_in);
            nmeas = length(specs.ctrl_out);
            P = simplify(P.content(1));
            [self.K,self.info] = LPV_unstructured_mix(P,nmeas,ncont,alpha,channels,P.parameters(),self.options);
            self.solved = true;
            
            % Save solver output
            self.performance = specs.performance;
            self.gamma = zeros(1,length(self.performance));
            self.mu = zeros(1,length(self.performance));
            
            for k = 1:length(channels)
                if k <= specs.nobj
                    val = sqrt(self.info.objective);
                else
                    val = 1;
                end
                
                if channels(k).performance == Inf
                    self.gamma(k) = val;
                else
                    self.mu(k) = val;
                end
            end
            
            if strcmp(self.options.objective,'wc')
                gamormu = self.gamma + self.mu;
                if specs.nobj > 0
                    obj = 1:specs.nobj;
%                    self.performance(obj) = self.performance(obj).*(1./gamormu(obj,1));
                end
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
            cap.parametric = true;
            cap.fixedorder = false;
        end
    end
end