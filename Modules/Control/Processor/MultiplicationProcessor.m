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

classdef MultiplicationProcessor
    %MULTIPLICATIONPROCESSOR Processing based on pre/post multiplication of
    %the controller
    %   The generalized plant is pre- and postmultiplied by Q, resp. R. Q 
    %   is used to make unstable modes controllable from the generalized 
    %   plant's perspective. R should render the plant proper when improper
    %   weights are used that make the generalized plant improper.
    %   The controller is retrieved by multiplying with Q and R:
    %   Kreal = Q*Kvirt*R
    
    properties
        Q;  % make the plant stabilizable
        R;  % make the plant proper
    end
    
    methods
        function [self,configuration2,problem2,solverhint] = preprocess(self,configuration,problem)
            
            % Build generalized plant from input data
            generalizedplant = constructGeneralizedPlant(problem,configuration);
            generalizedplant = build(generalizedplant);
            P = generalizedplant.sys_G;
            
            % Set some constants
            [nz,nw] = exogsize(P);
            [nmeas,ncont] = ctrlsize(P);
            P = ss(P);
            
            % SISO CASE
            P11 = P(1:nz,1:nw);
            P12 = P(1:nz,nw+(1:ncont));
            P21 = P(nz+(1:nmeas),1:nw);

            % Set weights
            self.Q = 1;
            self.R = 1;

            % Check for the degree of improperness (# - #)
            [z,p] = zpkdata(P12);
            improper = max(cellfun(@(x,y)(length(x)-length(y)),z(:),p(:)))
            if improper > 0
                p = pole(P);
                p = max(abs([p(p<eps^-0.5);0.1]))*10;
                self.R = ss(tf(1,butterpoly(improper,p)));
            end

            % Check for unstable weights (# - #)
            p11 = pole(minreal(P11)); % fix ill conditioning
            p21 = pole(minreal(P21));
            
            unstable = sum(abs(p11)<=sqrt(eps)) - sum(abs(p21)<=sqrt(eps))
            if unstable > 0
                den = zeros(1,unstable + 1);
                den(1) = 1;
                z = min(abs([p11(p11>sqrt(eps));p21(p21>sqrt(eps));10]))/10; % Maybe look at all poles in the system?
                self.Q = ss(tf(butterpoly(unstable,z),den));
            end
            
            % Adjust configuration
            configuration2 = configuration;
            Rl = LTIsys(self.R,Signal(size(self.R,2)),configuration2.ctrl_in);
            Ql = LTIsys(self.Q,configuration2.ctrl_out,Signal(size(self.Q,1)));
            configuration2 = configuration2.addsystem(Rl);
            configuration2 = configuration2.addsystem(Ql);
            configuration2.ctrl_in = Rl.in;
            configuration2.ctrl_out = Ql.out;
            
            % Adjust the problem: in and outputs
            problem2 = problem;
            problem2.ctrl_in = Rl.in;
            problem2.ctrl_out = Ql.out;
            
            % Add a solverhint so that the decision for the solver is
            % forced in some direction
            solverhint.unstable = false;
            solverhint.improper = false;
        end
        
        function [self,K] = postprocess(self,K)
            % Reconstruct the controller
            K = self.Q*K*self.R;
        end
    end
    
end

