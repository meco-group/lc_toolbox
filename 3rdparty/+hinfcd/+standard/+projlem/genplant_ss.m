classdef genplant_ss
% GENPLANT_SS A generalized plant with helper functions for the projection 
% lemma approach with order reductions in state-space form
%
% This class mainly provides convenient functions to obtain subparts of the
% state-space realization matrices that appear in the LMIs of the controller
% synthesis problem. As such, the implementation of the actual synthesis
% becomes easier to read and understand. 

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
        A      % the plant's A matrix
        Bu     % the plant's Bu matrix
        Bw     % the plant's Bw matrix
        Cy     % the plant's Cy matrix
        Cz     % the plant's Cz matrix
        Dyw    % the plant's Dyw matrix
        Dzu    % the plant's Dzu matrix
        Dzw    % the plant's Dzw matrix 
        Dyu1   % the plant's Dyu matrix
        Dyu2   % the plant's Dyu matrix after the elimination of the impulsive modes
        Duy    % the feedback gain for regularizing the descriptor system 
        r      % number of "direct" controller order reductions possible through unconstrained inputs
        s      % number of "direct" controller order reductions possible through measured outputs
        lambda % unstable or infinite invariant zero of the original plant that gives rise to order reductions
        i      % 0 if lambda is an invariant zero of Pzu, 1 if it is an invariant zero of Pyw
        Ts     % sampling time
    end
    
    methods
        function obj = genplant_ss(prob)
            % GENPLANT_SS Constructor  
            
            % load the plant
            Piobar = prob.Piobar();
            obj.Ts = Piobar.Ts; 
            nr = rank(Piobar.E);
            ny = prob.ny() + prob.no(); 
            nu = prob.nu() + prob.ni(); 
                       
            % transform to state-space
                % partition plant matrices
                A11 = Piobar.A(1:nr,1:nr);              A12 = Piobar.A(1:nr,(nr+1):end);            Bw1 = Piobar.B(1:nr,1:(end-nu));            Bu1 = Piobar.B(1:nr,(end-nu+1):end); 
                A21 = Piobar.A((nr+1):end,1:nr);        A22 = Piobar.A((nr+1):end,(nr+1):end);      Bw2 = Piobar.B((nr+1):end,1:(end-nu));      Bu2 = Piobar.B((nr+1):end,(end-nu+1):end);
                Cz1 = Piobar.C(1:(end-ny),1:nr);        Cz2 = Piobar.C(1:(end-ny),(nr+1):end);      Dzw = Piobar.D(1:(end-ny),1:(end-nu));      Dzu = Piobar.D(1:(end-ny),(end-nu+1):end);
                Cy1 = Piobar.C((end-ny+1):end,1:nr);    Cy2 = Piobar.C((end-ny+1):end,(nr+1):end);  Dyw = Piobar.D((end-ny+1):end,1:(end-nu));  obj.Dyu1 = Piobar.D((end-ny+1):end,(end-nu+1):end);
                
                % calculate regularizing static output feedback gain
                % See V. Lovass-Nagy, D.L. Powers, R.J. Schilling. "On
                % regularizing descriptor systems by output feedback.", IEEE
                % Transactions on Automatic Control, vol. 39, no. 7, 1994.
                ns = rank(A22,1e-8); 
                [U,S,V] = svd(A22); 
                S((ns+1):end,(ns+1):end) = eye(length(S)-ns); U = U*sqrt(S); V = sqrt(S)*V'; 
                F = U\Bu2; F2 = F((ns+1):end,:); 
                G = Cy2/V; G2 = G(:,(ns+1):end);
                obj.Duy = -pinv(F2)*pinv(G2);
            
                % explicit LFT calculation
                A11 = A11+Bu1*obj.Duy*Cy1;  A12 = A12+Bu1*obj.Duy*Cy2;  Bw1 = Bw1+Bu1*obj.Duy*Dyw;
                A21 = A21+Bu2*obj.Duy*Cy1;  A22 = A22+Bu2*obj.Duy*Cy2;  Bw2 = Bw2+Bu2*obj.Duy*Dyw; 
                Cz1 = Cz1+Dzu*obj.Duy*Cy1;  Cz2 = Cz2+Dzu*obj.Duy*Cy2;  Dzw = Dzw+Dzu*obj.Duy*Dyw; 
                
                % explicit algebraic/nondynamic mode elimination
                AA = A22\A21;
                ABu = A22\Bu2;
                ABw = A22\Bw2;
                
                obj.A = A11-A12*AA;     obj.Bw = Bw1-A12*ABw;   obj.Bu = Bu1-A12*ABu;
                obj.Cz = Cz1-Cz2*AA;    obj.Dzw = Dzw-Cz2*ABw;  obj.Dzu = Dzu-Cz2*ABu; 
                obj.Cy = Cy1-Cy2*AA;    obj.Dyw = Dyw-Cy2*ABw;  obj.Dyu2 = -Cy2*ABu; 
                
                % balance
                Piobar = balreal(ss(obj.A,[obj.Bw, obj.Bu],[obj.Cz ; obj.Cy],[obj.Dzw obj.Dzu ; obj.Dyw obj.Dyu2])); 
                % Piobar = (ss(obj.A,[obj.Bw, obj.Bu],[obj.Cz ; obj.Cy],[obj.Dzw obj.Dzu ; obj.Dyw obj.Dyu2])); 
                Piobar.E = eye(size(Piobar.A)); 
                obj.A = Piobar.A;                       obj.Bw = Piobar.B(:,1:(end-nu));                    obj.Bu = Piobar.B(:,(end-nu+1):end);
                obj.Cz = Piobar.C(1:(end-ny),:);        obj.Dzw = Piobar.D(1:(end-ny),1:(end-nu));          obj.Dzu = Piobar.D(1:(end-ny),(end-nu+1):end);
                obj.Cy = Piobar.C((end-ny+1):end,:);    obj.Dyw = Piobar.D((end-ny+1):end,1:(end-nu));      obj.Dyu2 = Piobar.D((end-ny+1):end,(end-nu+1):end);
                
            % determine number of possible order reductions
            Piobar.D((end-ny+1):end,(end-nu+1):end) = 0; 
            [~,nr,i,lambda] = hinfcd.util.orderbound(Piobar,prob.ny()+prob.no(),prob.nu()+prob.ni());
            if i == 0
                obj.r = nr;
                obj.s = rank(prob.Wo.E); 
            else
                obj.r = rank(prob.Wi.E); 
                obj.s = nr;
            end
            obj.lambda = lambda; 
            obj.i = i;

        end
        
        %% dimensions
        function ny = ny(obj)
        % NY Returns the number of measured outputs
            ny = size(obj.Cy,1);
        end
        
        function nu = nu(obj)
        % NU Returns the number of control inputs
            nu = size(obj.Bu,2);
        end
        
        function nz = nz(obj)
        % NW Returns the number of performance inputs
            nz = size(obj.Cz,1);
        end
        
        function nw = nw(obj)
        % NW Returns the number of performance inputs
            nw = size(obj.Bw,2);
        end
        
        function nx = nx(obj)
        % NX Returns the number of states
            nx = size(obj.A,1);
        end
        
        function sX = sizeX(obj)
        % SX Returns the dimensions of the Lyapunovmatrix X
            sX = [obj.nx() obj.nx()-obj.s];
        end
        
        function sY = sizeY(obj)
        % SY Returns the dimension sof the Lyapunovmatrix Y
            sY = [obj.nx()-obj.r obj.nx()];
        end
        
        %% matrices that are involved in the performance LMI with X
        function NX = NX(obj)
        % NX Returns the outer factor of the performance LMI with X
            U = null([obj.Cy(1:obj.ny(),1:obj.nx()) obj.Dyw(1:obj.ny(),:)]);
            U((obj.nx()+1):obj.nx(),:) = 0;
            NX = blkdiag(U,eye(obj.nz()));
        end

        %% matrices that are involved in the performance LMI with Y
        function NY = NY(obj)
        % NY Returns the outer factor of the performance LMI with Y
            V = null([obj.Bu([1:obj.nx(),(obj.nx()+1):obj.nx()], 1:obj.nu())' obj.Dzu(:,1:obj.nu())']);
            NY = blkdiag(V,eye(obj.nw()));
        end
        
    end
end

