function [obj,P,Q,part] = kalreal(obj,zerotol)
% KALREAL Returns the Kalman canonical form of a descriptor system
% Returns the Kalman canonical form of a descriptor system and gives the
% dimensions of the 4 different "subsystems". This function implements the
% geometric approach from A. Banaszuk, M. Kociecki and F.L. Lewis,
% "Kalman Decomposition for Implicit Linear Systems", IEEE Transactions on
% Automatic Control, vol. 37, no. 10, October 1992. 

% b = KALREAL(a) returns the canonical realization b corresponding the 
% the Kalman decomposition of a
%
% [b,u,v] = KALREAL(a) returns the same as the previous syntax, but
% additionally returns two orthogonal matrices u and v that transform
% realization a into its Kalman decomposition (u*A*v, u*B, C*v, D, u*E*v)
%
% [b,u,v,part] = KALREAL(a) returns the same as the previous syntax, and 
% the dimensions of the subsystems, e.g. part.nro is the number of states 
% that are not reachable (nr) but observable (o)
%
% [...] = KALREAL(a,zerotol) returns the same as the previous syntax, but
% additionaly considers numbers that are smaller than zerotol (default: 
% sqrt(eps)) in absolute value as zero while calculating a basis for the 
% required subspaces
%
% See also HINFCD.DSS.MINREAL. 

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

    import hinfcd.util.*;
    if nargin<2
        zerotol = sqrt(eps);
    else
        assert(isnumeric(zerotol) && isscalar(zerotol) && zerotol>=0, 'KALREAL: zerotol must be a positive number.'); 
    end
    
    % 0. Check for static systems
    if isempty(obj.E) && isempty(obj.A)
        P = []; Q = [];
        part = struct(); 
        part.rno = 0; part.ro = 0; part.nrno = 0; part.nro = 0;
        return;
    end

    % 1. Calculation of the intermediate geometric subspaces
    n = size(obj.E,1);
    VB = wongseq(obj.E,obj.A,obj.B,zeros(n),eye(n));
    RB = wongseq(obj.A,obj.E,obj.B,zeros(n),zeros(n));
    VC = wongseq(obj.E,obj.A,zeros(n),obj.C,eye(n));
    RC = wongseq(obj.A,obj.E,zeros(n),obj.C,zeros(n));
    
    % 2. Calculation of the observability and controllability subspaces
    CB = isectspace(RB,VB);     % CB = intersection of RB and VB 
    OC = orth([VC,RC]);         % OC = sum of VC and RC
    
    % 3. Construction of the transformation matrices
        % 3.1 First columns
        X1 = isectspace(CB,OC);                        % X1 = intersection of CB and OC
        Z1 = orthtol([obj.E*X1,obj.A*X1],zerotol);     % Z1 = sum of E(X1) and A(X1)
        
        % 3.2 Other columns
        X2 = dirsumdec(X1,CB);
        X3 = dirsumdec(X1,OC);
        X4 = dirsumdec(X1,X2,X3,eye(n));
        
        Z2 = dirsumdec(Z1,orthtol([obj.E*CB,orth(obj.B)],zerotol));
        Z3 = dirsumdec(Z1,orthtol([obj.E*OC,obj.A*OC],zerotol));
        Z4 = dirsumdec(Z1,Z2,Z3,eye(n));
        
        % 3.3 Build P and Q
        P = [Z1 Z2 Z3 Z4]';
        Q = [X1 X2 X3 X4];
        
    % 4. Return subparts of the system
    part = struct(); 
    part.rno = size(X1,2);  % reachable but not observable
    part.ro = size(X2,2);   % both reachable and observable 
    part.nrno = size(X3,2); % not reachable nor observable
    part.nro = size(X4,2);  % not reachable but observable

    % 5. Transform
    if any(part.rno) || any(part.nro) || any(part.nrno)
        Ek = P*obj.E*Q; Ek = roundsmall(Ek,zerotol);
        Ak = P*obj.A*Q; Ak = roundsmall(Ak,zerotol); 
        Bk = P*obj.B;   Bk = roundsmall(Bk,zerotol);
        Ck = obj.C*Q;   Ck = roundsmall(Ck,zerotol); 
        obj.set('E',Ek,'A',Ak,'B',Bk,'C',Ck);
    else
        % do not do the transformation to save numerical troubles
        P = eye(size(obj.E)); Q = eye(size(obj.E));
    end

end
