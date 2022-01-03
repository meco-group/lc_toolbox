function obj = minreal(obj,zerotol)
% MINREAL Minimal realization of a descriptor system
%
% If DSTOOLS is found on the path, its implementation is used; 
% see A. Varga, "DSTOOLS - The Descriptor System Tools for MATLAB.", 2019.
% If not, the minimal realization is found in two steps:
%   1. A Kalman decomposition, eliminating the unreachable and unobservable
%   states, resulting in a 'conditionally' minimal realization. 
%   2. A Kronecker-Weierstrass decomposition, bringing the 'conditionally'
%   minimal realization into the 'deflated' minimal realization.
% This procedure is based on chapter 4 of V.I. Sokolov, P. Benner (sup.), 
% A.C. Antoulas (sup.) and V. Mehrmann (sup)., "Contributions to the 
% Minimal Realization Problem for Descriptor Systems", PhD thesis, Faculty
% of Mathematics, TU Chemnitz, Germany, 2006. 
%
% Note: MINREAL removes uncontrollable and unobservable states without
% notice. Removing uncontrollable states may effect the output.
%
% b = MINREAL(a) returns the minimal realization b of a
%
% b = MINREAL(a,zerotol) returns the same as the previous syntax, but
% additionaly considers numbers that are smaller than zerotol in absolute 
% value as zero while eliminating unobservable, uncontrollable or 
% redundant states
%
% See also HINFCD.DSS.KALREAL, HINFCD.DSS.KRONREAL. 

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

    if nargin<2
        zerotol = sqrt(eps);
    else
        assert(isnumeric(zerotol) && isscalar(zerotol) && zerotol>=0, 'MINREAL: zerotol must be a positive number.'); 
    end
     
    if exist('gminreal','file') % DSTOOLS implementation
        if nargin<2; mr = gminreal(obj,0); else; mr = gminreal(obj,zerotol); end
        if isempty(mr.E); mr.E = eye(size(mr.A)); end
        obj.set('A',mr.A,'B',mr.B,'C',mr.C,'D',mr.D,'E',mr.E);
        
    else % own naive implementation
        % 1. Kalman decomposition
        %    -> removes unobservable and uncontrollable states 
        %       (cf. state-space realizations)
        [obj,~,~,part] = kalreal(obj,zerotol);
        Ered = obj.E(part.rno+1:part.rno+part.ro,part.rno+1:part.rno+part.ro);
        Ared = obj.A(part.rno+1:part.rno+part.ro,part.rno+1:part.rno+part.ro);
        Bred = obj.B(part.rno+1:part.rno+part.ro,:); 
        Cred = obj.C(:,part.rno+1:part.rno+part.ro); 
        obj = obj.set('E',Ered,'A',Ared,'B',Bred,'C',Cred);

        % 2. Kronecker-Weierstrass decomposition
        %    -> removes redundant algebraic states
        [obj,~,~,feigv] = kronreal(obj,zerotol);
        [Js,V,nidxone] = hinfcd.util.niljordan(obj.E((feigv+1):end,(feigv+1):end));
        E = obj.E; B = obj.B; C = obj.C; D = obj.D;
        E((feigv+1):end,(feigv+1):end) = Js;
        B((feigv+1):end,:) = V\obj.B((feigv+1):end,:);
        C(:,(feigv+1):end) = obj.C(:,(feigv+1):end)*V;
        E = E(1:end-nidxone,1:end-nidxone); 
        D = D-C(:,(feigv+1):end)*B((feigv+1):end,:);
        B = B(1:end-nidxone,:);
        C = C(:,1:end-nidxone); 
        A = obj.A(1:end-nidxone,1:end-nidxone);
        obj.set('E',E,'A',A,'B',B,'C',C,'D',D);

        % display results
        fprintf('MINREAL: Removed %i unobservable and/or uncontrollable state(s).\n', part.rno+part.nrno+part.nro);
        fprintf('MINREAL: Removed %i redundant algebraic state relation(s).\n', nidxone);
    end

end
