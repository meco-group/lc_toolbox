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

classdef Solver_systune < Solver
    %SOLVER_SYSTUNE Summary of this class goes here
    %   Detailed explanation goes here
        
    methods
        function self = Solver_systune(~)
            % Do nothing for now..
        end
        
        function self = solve(self,config,specs,vars)
            % Compute plant state-space
            if ~(all(scale(specs.performance(1:specs.nobj)))==1)
                warning('Systune rescales the input weights by the provided scaling factors');
            end
            specs = specs.rescale('obj');
            [P,wspecs] = Solver.plant(config,specs,vars);
            ch = Solver.channels(wspecs);   

            % Split up the system in parts
            ncont = length(specs.ctrl_in);
            nmeas = length(specs.ctrl_out);
            nw = 1:(size(P,2)-ncont);
            nu = (nw(end)+1):size(P,2);
            nz = 1:(size(P,1)-nmeas);
            ny = (nz(end)+1):size(P,1);

            % define systune tuning goals
            if specs.nobj>1
                warning('Systune cannot handle a true multi-objective problem. For multiple objectives, it minimizes max(norm(T_i)) instead of sum(norm(T_i)).');
            end            
            P = std(P);
            P.inputname = arrayfun(@(x)['u' num2str(x)],1:size(P,2),'UniformOutput',false);
            P.outputname = arrayfun(@(x)['y' num2str(x)],1:size(P,1),'UniformOutput',false);
            Req = cell(1,length(ch.In));
            for k = 1:length(ch.In)
                inputs = arrayfun(@(x)['u' num2str(x)],find(sum(ch.In{k},2)),'UniformOutput',false);
                outputs = arrayfun(@(x)['y' num2str(x)],find(sum(ch.Out{k},1)),'UniformOutput',false);
                
                sc = 1/scale(specs.performance(k));
                if ismember(k,ch.H2)
                    Req{k} = TuningGoal.Variance(inputs,outputs,sc);
                else %Hinf
                    Req{k} = TuningGoal.Gain(inputs,outputs,sc);
                end
            end
                        
            objectives = horzcat(Req{1:specs.nobj});
            constraints = horzcat(Req{specs.nobj+1:end});

            % make controller block
            if specs.order == -1
                K = ltiblock.ss('K',size(P.a,1),length(ny),length(nu));
            else
                K = ltiblock.ss('K',specs.order,length(ny),length(nu));
            end
            K.u = P.outputname(ny);
            K.y = P.inputname(nu);

            % Try solving for a controller
            P = connect(P,K,P.inputname(nw),P.outputname(nz));
            [CL,fSoft,gHard,info] = systune(P,objectives,constraints);
            gammaormu = [fSoft,gHard];
            
            % rescale performance
            self.performance = specs.performance;
            if specs.nobj > 0
                obj = 1:specs.nobj;
                self.performance(obj) = dealscale(self.performance(obj),1./fSoft(obj));
            end
            
            % save output
            self.K = fromstd(getBlockValue(CL,'K'));
            self.gamma = gammaormu;
            self.gamma(ch.H2) = 0;
            self.mu = gammaormu;
            self.mu(ch.Hinf) = 0;
            self.info = mergestruct(self.info,info);
            self.solved = true;
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

