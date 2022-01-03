classdef(InferiorClasses = {?ss}) dss < ss
% DSS Extension of standard MATLAB dss models
%
% See also DSS, SS.

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
    
    methods
        function obj = dss(varargin)
        % DSS Construct a state-space model using the standard MATLAB function
            if nargin<5
                error('Not enough input arguments.');
            elseif nargin==5
                v = varargin(1:4);
            else
                v = varargin([1:4,6:end]);
            end
            obj@ss(v{:});
            obj.E = varargin{5};
            if isempty(obj.E)
                obj.E = eye(size(obj.A));
            end
        end
        
        obj = minreal(obj,zerotol);
        [obj,u,v,feigv] = kronreal(obj,zerotol);
        [obj,u,v,part] = kalreal(obj,zerotol);
        [obj,g,h] = balreal(obj);
        [obj,r,s,M,N,P,Q] = oredreal(obj,ny,nu);
    end
    
    methods(Access=public)
        function b = isregular(obj)
        % ISREGULAR Returns whether the descriptor system is regular or not
        
        % See theorem 1 of T. Berger et al., "New lower bound for the
        % distance to singularity of regular matrix pencils", Proceedings 
        % in Applied Mathematics and Mechanics, vol. 17, 2017.
            n = size(obj.E,1);
            E = repmat({obj.E},[n 1]);
            A = repmat({obj.A},[n 1]);
            W1 = blkdiag(A{:}); W1 = [W1 ; zeros(n,n^2)];
            W2 = blkdiag(E{:}); W2 = [zeros(n,n^2) ; W2];  
            W = W1+W2;
            b = (rank(W)==size(W,2)); 
        end
        
        function b = isfdobsv(obj)
        % ISFDOBSV Returns whether the descriptor system is finite dynamics controllable or not
        
        % See J. Sjoberg, "Descriptor Systems and Control Theory", Technical report: 
        % LiTH-ISY-R-2688, Linköping University, 2005.
            GIF = impsep(obj); 
            [~,~,~,part] = kalreal(GIF); 
            b = (part.ro+part.nro)==size(GIF.E,1);
        end
        
        function b = isfddetect(obj)
        % ISFDDETECT Returns whether the descriptor system is finite dynamics detectable or not
        
        % See J. Sjoberg, "Descriptor Systems and Control Theory", Technical report: 
        % LiTH-ISY-R-2688, Linköping University, 2005.
            GIF = impsep(obj);
            [GIF,~,~,part] = kalreal(GIF);
            E = GIF.E([1:part.rno,(part.rno+part.ro+1):(part.rno+part.ro+part.nrno)],[1:part.rno,(part.rno+part.ro+1):(part.rno+part.ro+part.nrno)]);
            A = GIF.A([1:part.rno,(part.rno+part.ro+1):(part.rno+part.ro+part.nrno)],[1:part.rno,(part.rno+part.ro+1):(part.rno+part.ro+part.nrno)]);
            b = all(eig(A,E)<0);
        end
        
        function b = isfdctrb(obj)
        % ISFDCTRB Returns whether the descriptor system is finite dynamics controllable or not
        
        % See J. Sjoberg, "Descriptor Systems and Control Theory", Technical report: 
        % LiTH-ISY-R-2688, Linköping University, 2005.
            GIF = impsep(obj); 
            [~,~,~,part] = kalreal(GIF); 
            b = (part.ro+part.rno)==size(GIF.E,1);
        end
        
        function b = isfdstab(obj)
        % ISFDSTAB Returns whether the descriptor system is finite dynamics stabilizable or not
        
        % See J. Sjoberg, "Descriptor Systems and Control Theory", Technical report: 
        % LiTH-ISY-R-2688, Linköping University, 2005.
            GIF = impsep(obj);
            [GIF,~,~,part] = kalreal(GIF);
            E = GIF.E((part.rno+part.ro+1):end,(part.rno+part.ro+1):end);
            A = GIF.A((part.rno+part.ro+1):end,(part.rno+part.ro+1):end);
            b = all(eig(A,E)<0);
        end
        
        function b = isimpobsv(obj)
        % ISIMPOBSV Returns whether the descriptor system is finite dynamics observable or not
        % Finite dynamics observability is also referred to as 
        % R-observability.
        
        % See J. Sjoberg, "Descriptor Systems and Control Theory", 
        % Technical report: LiTH-ISY-R-2688, Linköping University, 2005.
            n = size(obj.E,1);
            b = rank([obj.E obj.A ; zeros(n) obj.E ; zeros(size(obj.C,1),n) obj.C])==(n+rank(obj.E));
        end
        
        function b = isimpctrb(obj)
        % ISIMPCTRB Returns whether the descriptor system is impulse controllable or not
        
        % See J. Sjoberg, "Descriptor Systems and Control Theory", 
        % Technical report: LiTH-ISY-R-2688, Linköping University, 2005.
            n = size(obj.E,1);
            b = rank([obj.E zeros(n) zeros(n,size(obj.B,2)) ; obj.A obj.E obj.B])==(n+rank(obj.E));  
        end
        
        function [GS,GNS] = stabsep(varargin)
        % STABSEP Separates the stable and the unstable dynamics
        % See also STABSEP. 
            assert(isa(varargin{1},'hinfcd.dss'), 'First argument of stabsep should be a hinfcd.dss object.'); 
            [obj,~,~,nslow] = kronreal(varargin{1}); 
            Gslow = hinfcd.dss(obj.A(1:nslow,1:nslow),obj.B(1:nslow,:),obj.C(:,1:nslow),obj.D,obj.E(1:nslow,1:nslow),obj.Ts);
            Gfast = hinfcd.dss(obj.A(nslow+1:end,nslow+1:end),obj.B(nslow+1:end,:),obj.C(:,nslow+1:end),zeros(size(obj.D)),obj.E(nslow+1:end,nslow+1:end),obj.Ts);
            [GS,GNS] = stabsep@DynamicSystem(Gslow,varargin{2:end}); 
            GS = GS + Gfast; 
            GS = hinfcd.dss(GS.A,GS.B,GS.C,GS.D,GS.E,GS.Ts);
            GNS = hinfcd.dss(GNS.A,GNS.B,GNS.C,GNS.D,GNS.E,GNS.Ts);
        end
        
        function [GIF,GIMP] = impsep(obj)
        % IMPSEP Separates the finite dynamics and the impulsive dynamics
        % In fact, this is nothing but separating the slow and the fast
        % subsystem. 
            [obj,~,~,nslow] = kronreal(obj);
            GIF = hinfcd.dss(obj.A(1:nslow,1:nslow),obj.B(1:nslow,:),obj.C(:,1:nslow),obj.D,obj.E(1:nslow,1:nslow),obj.Ts);
            GIMP = hinfcd.dss(obj.A(nslow+1:end,nslow+1:end),obj.B(nslow+1:end,:),obj.C(:,nslow+1:end),zeros(size(obj.D)),obj.E(nslow+1:end,nslow+1:end),obj.Ts);
        end
        
    end

end
    