function [obj,r,s,T,P] = oredrealsvd(obj,ny,nu,zerotol)
% OREDREALSVD SVD canonical realization suited for reduced-order controller synthesis
% 
% b = OREDREALSVD(a,ny,nu) returns b, an equivalent realization of a, with
% one of the following canonical forms:
%
% 1. if r > s: 
%     [Ie 0]*dx     =   Ab*xb  + Bwb*w  + [Bub1   0]*ub
%     [0  0]                              [Bub2  Ir]
%                                         [Bub3   0]
%             z     =   Czb*xb + Dzwb*w + [Dzub1  0]*ub
%             yb    =   Cyb*xb + Dywb*w + Dyub*ub
%  
%     In this case, s is set to 0.
%  
% 2. if s >= r:
%     [Ie 0]*dx     =   Ab*xb               + Bwb*w     + Bub*ub
%     [0  0] 
%             z     =   Czb*xb              + Dzwb*w    + Dzub*ub
%             yb    =   [Cyb1 Cyb2 Cyb3]*xb + [Dywb1]*w + Dyub*ub
%                       [0      Is    0]      [    0]
%
%     In this case, r is set to 0.
%
% ny is the dimension of the measured outputs y, nu is the dimension of the
% control inputs u. r is the achieved size of Ir (case 1), whereas s is the 
% achieved size of Is (case 2). Ie is the unity matrix of size e = rank(E). 
% Ir and Is are identity matrices. Bub3 has as n-e rows (case 1), whereas 
% Cyb3 has n-e columns (case 2). 
% IMPORTANT: The new plant has scaled control inputs ub or scaled measured
% outputs yb, i.e. the original signals are ua = P*ub or ya = P*yb 
% (see next syntax). 
%
% [b,r,s,T,P] = OREDREALSVD(a,ny,nu) returns the same as the previous 
% syntax, and furthermore returns the transformation matrices that were 
% used to obtain the canonical form as follows: xa = T*xb and ua = P*ub in
% case 1 or ya = P*yb in case 2, where 'a' and 'b' refer to 
% the original and the returned realization respectively. 
%
% See also HINFCD.UTIL.OREDREAL. 

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
    if nargin<4
        zerotol = sqrt(eps);
    else
        assert(isnumeric(zerotol) && isscalar(zerotol) && zerotol>=0, 'OREDREAL: zerotol must be a positive number.'); 
    end
    
    % 1. Calculate the SVD canonical form
    re = rank(obj.E);
    [U,S,V] = svd(obj.E);
    U = U*blkdiag(sqrt(S(1:re,1:re)),eye(size(obj.E)-re));
    V = V*blkdiag(sqrt(S(1:re,1:re)),eye(size(obj.E)-re));
    A = U\obj.A/V'; 
    B = U\obj.B;
    C = obj.C/V';
    obj = obj.set('A',A,'B',B,'C',C,'E',diag([ones(1,re) zeros(1,size(obj.E,1)-re)]));
    
    % 2. Partition the plant
    assert(nu <= size(obj.B,2), 'OREDREAL: nu larger than number of inputs.');
    assert(ny <= size(obj.C,1), 'OREDREAL: ny larger than number of outputs.');
    E = obj.E;      A = obj.A;                      Bw = obj.B(:,1:(end-nu));               Bu = obj.B(:,(end-nu+1):end);
                    Cz = obj.C(1:(end-ny),:);       Dzw = obj.D(1:(end-ny),1:(end-nu));     Dzu = obj.D(1:(end-ny),(end-nu+1):end);
                    Cy = obj.C((end-ny+1):end,:);   Dyw = obj.D((end-ny+1):end,1:(end-nu)); Dyu = obj.D((end-ny+1):end,(end-nu+1):end);
    n = size(obj.A,1);

    % 3. Determine all possible order reductions 
    rr = rank(Bu(1:re,:)); r = min([nu-rank(Dzu),rr]);
    ss = rank(Cy(:,1:re)); s = min([ny-rank(Dyw),ss]); 

    % 4. Calculate transformation matrices
         if r>s
         % Case 1: r > s
            r = min([r,re]); s = 0;
            
            % input scaling
            [~,~,Q] = svd(Dzu);

            % state transformation
            Bub = Bu*Q;
            [~,~,M] = svd(Bub');
            M = [M(:,(rr+1):end) M(:,1:rr)]; BubM = Bub'*M;
            M(:,(end-rr+1):end) = M(:,(end-rr+1):end)/BubM((end-rr+1):end,(end-rr+1):end);
            M = circshift(M',re,1);
            T = inv(M);
            P = eye(size(obj,1));
            Q = blkdiag(eye(size(obj,2)-nu),Q);
            
         else
         % Case 2: s >= r
            r = 0; s = min([s,re]);
            
            % input scaling
            [~,~,M] = svd(Dyw'); P = M'; 

            % state transformation
            Cyb = P*Cy(:,1:re); 
            [~,~,N] = svd(Cyb);
            N = [N(:,(ss+1):end) N(:,1:ss)]; CybN = Cyb*N;
            N(:,(end-ss+1):end) = N(:,(end-ss+1):end)/CybN((end-ss+1):end,(end-ss+1):end);
            T = blkdiag(N,eye(n-re));
            P = blkdiag(eye(size(obj,1)-ny),P);
            Q = eye(size(obj,2)); 
            
         end

    % 4. Transform remaining matrices and set new realization
    Ab = T\obj.A*T; % pseudoinverse!
    Bb = T\obj.B*Q;
    Cb = P*obj.C*T;
    Db = P*obj.D*Q;

    Ab = roundsmall(Ab,zerotol);
    Bb = roundsmall(Bb,zerotol);
    Cb = roundsmall(Cb,zerotol);
    Db = roundsmall(Db,zerotol);

    obj = obj.set('A',Ab,'B',Bb,'C',Cb,'D',Db);
    if r>s; P = Q; end

end
