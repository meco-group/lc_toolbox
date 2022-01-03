classdef genplant_d
% GENPLANT A generalized plant with helper functions for the projection 
% lemma approach with order reductions in descriptor form
%
% This class mainly provides convenient functions to obtain subparts of the
% descriptor realization matrices that appear in the LMIs of the controller
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
        E      % the plant's E matrix
        r      % number of "direct" controller order reductions possible through unconstrained inputs
        s      % number of "direct" controller order reductions possible through measured outputs
        lambda % unstable or infinite invariant zero of the original plant that gives rise to order reductions
        i      % 0 if lambda is an invariant zero of Pzu, 1 if it is an invariant zero of Pyw
        rmax   % maximal number of controller orders that can be eliminated from the performance LMIs
        smax   % maximal number of controller orders that can be eliminated from the performance LMIs
        Ts     % sampling time
    end
    
    methods
        function obj = genplant_d(prob,real)
            % GENPLANT_D Constructor  
            
            % determine number of possible order reductions
            [~,nr,i,lambda] = hinfcd.util.orderbound(prob.P,prob.ny(),prob.nu());
            obj.s = rank(prob.Wo.E)+(i==1)*nr; 
            obj.r = rank(prob.Wi.E)+(i~=1)*nr;
            obj.lambda = lambda;
            obj.i = i;
            
            if nargin>1 && strcmp(real,'bal')
                % balanced realization in SVD form -> numerical attractive form for reconstruction
                Piobar = svdreal(balreal(prob.Piobar())); 
                rmax = 0;
                smax = 0;
            else     
                % canonical realization for order reductions -> numerical attractive form for synthesis
                [Piobar,rmax,smax] = oredreal(prob.Piobar(),prob.ny()+prob.no(),prob.nu()+prob.ni());
                % Piobar = svdreal(balreal(prob.Piobar())); 
                % rmax = 0;
                % smax = 0;
            end

            % variables that could be explicitly eliminated from the performance LMIs
            obj.smax = smax; 
            obj.rmax = rmax; 

            % partition the plant
            ny = prob.ny() + prob.no(); 
            nu = prob.nu() + prob.ni(); 

            obj.E = Piobar.E; 
            obj.A = Piobar.A;                       obj.Bw = Piobar.B(:,1:(end-nu));                    obj.Bu = Piobar.B(:,(end-nu+1):end);
            obj.Cz = Piobar.C(1:(end-ny),:);        obj.Dzw = Piobar.D(1:(end-ny),1:(end-nu));          obj.Dzu = Piobar.D(1:(end-ny),(end-nu+1):end);     
            obj.Cy = Piobar.C((end-ny+1):end,:);    obj.Dyw = Piobar.D((end-ny+1):end,1:(end-nu)); 
            
            obj.Ts = Piobar.Ts; 

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
        function AX = AX(obj)
        % AX Returns the submatrix of A for multiplication with X
        %
        % AX = [A11 ... A13]   if  r >= s
        %      [... ... ...]
        %      [A41 ... A43]
        %
        % AX = [A11 A12]   if r < s
        %      [... ...]
        %      [A41 A42]
        %
        % Note: AX = [AX1 ; AX2]
            AX = obj.A(:,1:(obj.nx()-obj.smax));
        end
        
        function CzX = CzX(obj)
        % CZX Returns the part of Cz that appears in the performance LMI with X
        %
        % CzX = [Cz1 ... Cz3]   if  r >= s
        %
        % CzX = [Cz1 Cz2]   if r < s
            CzX = obj.Cz(:,1:(obj.nx()-obj.smax)); 
        end
        
        function NX = NX(obj)
        % NX Returns the outer factor of the performance LMI with X
            U = null([obj.Cy(1:(obj.ny()-obj.smax),1:(obj.nx()-obj.smax)) obj.Dyw(1:(obj.ny()-obj.smax),:)]);
            U((obj.nx()-obj.smax+1):(obj.nx()-obj.smax),:) = 0;
            NX = blkdiag(U,eye(obj.nz()));
        end

        %% matrices that are involved in the performance LMI with Y
        function AY = AY(obj)
        % AY Returns the submatrix of A for multiplication with Y
        %
        % AY = [A11 ... A14]   if  r >= s
        %      [A41 ... A44]
        %
        % AY = [A11 ... A14]   if r < s
        %      [A21 ... A24]
        %      [A41 ... A44]
            AY = obj.A([1:(obj.nx()-obj.smax-obj.rmax) (obj.nx()-obj.smax+1):obj.nx()],:); 
        end
        
        function BwY = BwY(obj)
        % BWY Returns the part of Bw that appears in the performance LMI with Y
        %
        % BwY = [Bw1]   if r >= s
        %       [Bw4]
        %
        % BwY = [Bw1]   if r < s
        %       [Bw3]
        %       [Bw4]
            BwY = obj.Bw([1:(obj.nx()-obj.smax-obj.rmax), (obj.nx()-obj.smax+1):obj.nx()],:); 
        end
        
        function NY = NY(obj)
        % NY Returns the outer factor of the performance LMI with Y
            V = null([obj.Bu([1:(obj.nx()-obj.rmax-obj.smax),(obj.nx()-obj.smax+1):obj.nx()], 1:(obj.nu()-obj.rmax))' obj.Dzu(:,1:(obj.nu()-obj.rmax))']);
            NY = blkdiag(V,eye(obj.nw()));
        end
        
        %% matrices that are involved in the stability LMI
        function E1 = E1(obj)
        % E1 Returns the first part of the E-matrix involved in the stability LMI
        %
        % E1 = [E11 ... E13]    if r >= s
        %      [... ... ...]
        %      [E41 ... E43]
        %
        % E1 = [E11 E12]    if s < r
        %      [... ...]
        %      [E41 E42]
            E1 = obj.E(:,1:(obj.nx()-obj.s)); 
        end
        
        function E2 = E2(obj)
        % E2 Returns the second part of the E-matrix involved in the stability LMI
        %
        % E2 = [E11 ... E14]    if r >= s
        %      [E41 ... E44]
        %
        % E2 = [E11 ... E14]    if s < r
        %      [E31 ... E34]
        %      [E41 ... E44]
            E2 = obj.E([1:(obj.nx()-obj.r-obj.s), (obj.nx()-obj.s+1):obj.nx()],:); 
        end
        
        function I1 = I1(obj)
        % I1 Returns the first constant term of the stability LMI
        %
        % I1 is the top right constant term of the stability LMI, i.e.
        %   [E1' 0 ] [X   I1] >= 0 
        %   [0   E2] [I2  Y ]
            I1 = [blkdiag(eye(obj.nx()-obj.r-obj.s),zeros(obj.r,obj.s)) ; zeros(obj.s,obj.nx()-obj.r-obj.s) eye(obj.s)];
        end
        
        function I2 = I2(obj)
        % I2 Returns the second constant term of the stability LMI
        %
        % I2 is the bottom left constant term of the stability LMI, i.e.
        %   [E1' 0 ] [X   I1] >= 0 
        %   [0   E2] [I2  Y ]
        	I2 = [eye(obj.nx()-obj.s) ; zeros(obj.s, obj.nx()-obj.s)];
        end
        
        function NE12 = NE12(obj)
        % NE12 Returns the outer factor for the stability LMI to make it strict
        %
        % In case E1 and/or E2 are rank deficient, the LMI is not strict.
        % Since numerical solvers (esp. the ones based on interior-point
        % methods) have troubles with this, we can explicitly remove the
        % zero eigenvalues by an appropriate nonsquare outer factor NE12,
        % i.e.:
        %
        %   [E1' 0 ] [X   I1] >= 0   <=>   NE'*[E1' 0 ][X  I1]*NE >= 0
        %   [0   E2] [I2  Y ]                  [0   E2][I2 Y ]
           [U,~,~] = svd(obj.E1()'); M = U(:,1:rank(obj.E1()))';
           [U,~,~] = svd(obj.E2()); N = U(:,1:rank(obj.E2()))';
           NE12 = blkdiag(M,N)'; 
        end
        
        function [NE,NE1,NE2] = NE(obj)
        % NE Returns the outer factor for the stability LMI to make it strict
        %
        % In case E is rank deficient, the LMI is not strict.
        % Since numerical solvers (esp. the ones based on interior-point
        % methods) have troubles with this, we can explicitly remove the
        % zero eigenvalues by an appropriate nonsquare outer factor 
        % NE = blkdiag(NE1,NE2), i.e.:
        %
        %   [E' 0] [X  1] >= 0   <=>   NE'*[E' 0][X I]*NE >= 0
        %   [0  E] [I  Y]                  [0  E][I Y]
           [U,~,~] = svd(obj.E'); NE1 = U(:,1:rank(obj.E));
           [U,~,~] = svd(obj.E); NE2 = U(:,1:rank(obj.E));
           NE = blkdiag(NE1,NE2);
        end
        
    end
end

