classdef postprocessor 
% POSTPROCESSOR Postprocessing steps to the controller matrices originating from the standard problem

% This file is part of hinfcd.
% Copyright (c) 2019, Laurens Jacobs, MECO Research Team @ KU Leuven. 
% 
% hinfcd is free software: you can redistribute it and/or modify it under 
% the terms of the GNU Lesser General Public License as published by the 
% Free Software Foundation, version 3.
% 
% hinfcd is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for
% more details.
% 
% You should have received a copy of the GNU Lesser General Public License
% along with hinfcd. If not, see <https://www.gnu.org/licenses/>.
    
    properties
        % prob - The hinfcd problem
        % The extended H-infinity controller synthesis problem.
        prob
        
        % K - The controller resulting from the standard problem
        K
        
        % perf - Achieved performance
        perf
    end
    
    methods
        %% Constructor
        function obj = postprocessor(problem, K, perf)
        % POSTPROCESSOR Constructor
            
            % construct the object            
            assert(isa(problem,'hinfcd.problem'), 'Invalid problem.'); 
            obj.prob = problem; 
            
            assert(isa(K,'hinfcd.dss'), 'K must be a hinfcd.dss model.');
            obj.K = K; 
            
            obj.perf = perf;
            
        end
        
        %% Recombination the controller
        function obj = recombine(obj)
        % RECOMBINE Adds the impulsive and unstable modes of the weighting filters to the controller
            
            % controller parameters from the standard problem
            Ec = obj.K.E;
            Ac = obj.K.A; 
            Bc1 = obj.K.B(:,1:obj.prob.no()); Bc2 = obj.K.B(:,(obj.prob.no()+1):end);
            Cc1 = obj.K.C(1:obj.prob.ni(),:); Cc2 = obj.K.C((obj.prob.ni()+1):end,:); Cc = [Cc1 ; Cc2];
            Dc11 = obj.K.D(1:obj.prob.ni(),1:obj.prob.no()); Dc12 = obj.K.D(1:obj.prob.ni(),(obj.prob.no()+1):end);
            Dc21 = obj.K.D((obj.prob.ni()+1):end,1:obj.prob.no()); Dc22 = obj.K.D((obj.prob.ni()+1):end,(obj.prob.no()+1):end);
            Dc1 = [Dc11 ; Dc21]; Dc2 = [Dc12 ; Dc22];
            
            % separation variables
            [GAMMAo,~,Zotilde] = ofseparator(obj.prob); 
            [GAMMAi,~,Zi] = ifseparator(obj.prob);
            
            % reconstruct controller parameters
            E = blkdiag(obj.prob.Eo(),Ec,obj.prob.Ei()); 
            A = [obj.prob.Ao()-Zotilde*Dc1      -Zotilde*Cc     -(GAMMAo-Zotilde*Dc2)*Zi;
                 Bc1                            Ac              -Bc2*Zi                 ;
                 Dc11                           Cc1             obj.prob.Ai()-Dc12*Zi   ];
            B = [GAMMAo-Zotilde*Dc2;
                 Bc2               ;
                 Dc12              ];
            C = [Dc21 Cc2 GAMMAi-Dc22*Zi];
            D = Dc22;
            
            % return 
            obj.K = hinfcd.dss(A,B,C,D,E,obj.prob.Ts); 
        end
        
        %% Elimination of the nondynamic modes from the weighting filters
        function obj = elimstaticmodes(obj)
        % ELIMSTATICMODES Eliminates the static modes from the controller and returns a state-space model
        
            K = svdreal(obj.K); 
            nr = rank(K.E); 
            
            % partition matrices
            A11 = K.A(1:nr,1:nr);              A12 = K.A(1:nr,(nr+1):end);            B1 = K.B(1:nr,:); 
            A21 = K.A((nr+1):end,1:nr);        A22 = K.A((nr+1):end,(nr+1):end);      B2 = K.B((nr+1):end,:); 
            C1 = K.C(:,1:nr);                  C2 = K.C(:,(nr+1):end);                D = K.D;

            % check singularity of A22
            obj.prob.watchdog.warnCond(rank(A22)==length(A22),'Nondynamic modes of the controller are badly scaled. Results may be inaccurate.');
            
            % eliminate
            AA = A22\A21;
            AB = A22\B2;
            obj.K = set(obj.K,'A',A11-A12*AA,'B',B1-A12*AB,'C',C1-C2*AA,'D',D-C2*AB,'E',eye(nr));
            
        end
        
        %% Compensation for the feedthrough component of the generalized plant
        function obj = compensateDyu(obj)
        % COMPENSATEDYU Compensates for direct feedthrough in the generalized plant, if present
            
            % get synthesized controller
            [A0,B0,C0,D0] = dssdata(obj.K); 
            
            % get feedthrough matrix
            Dyu = obj.prob.Dyu(); 
            
            % compensate for the algebraic loop
            Mcyu = eye(size(obj.K,1))+D0*Dyu;
            C = Mcyu\C0;
            D = Mcyu\D0;
            A = A0-B0*Dyu*C;
            B = B0*(eye(size(obj.K,2))-Dyu*D); 
            
            % set updated controller data
            obj.K = set(obj.K,'A',A,'B',B,'C',C,'D',D);
            
        end
           
    end
    
end