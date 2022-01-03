function [obj,r,s,M,N,P,Q] = oredreal(obj,ny,nu,zerotol)
% OREDREAL Canonical realization suited for reduced-order controller synthesis
% 
% b = OREDREAL(a,ny,nu) returns b, an equivalent realization of a, with the 
% following canonical form:
%
%   Eb*dx =          Ab*x +     Bwb*w + [Bub1   0]*ub
%                                       [Bub2  Ir]
%                                       [Bub3   0]
%      z  =         Czb*x +    Dzwb*w + [Dzub1  0]*ub
%      yb = [Cyb1 Cyb2]*x + [Dywb1]*w +     [Dyub]*ub
%           [0      Is]     [    0]
%
% ny is the dimension of the measured outputs y, nu is the dimension of the
% control inputs u. r is the achieved size of Ir, whereas s is the achieved
% size of Is. Ir and Is are identity matrices. Bub3 has as s rows.
% IMPORTANT: The new plant has scaled control inputs ub and scaled measured
% outputs yb, i.e. the original inputs are ua = Q*ub and ya = P*yb 
% (see next syntax). 
%
% [b,r,s,M,N,P,Q] = OREDREAL(a,ny,nu) returns the same as the previous 
% syntax, and furthermore returns the transformation matrices that were 
% used to obtain the canonical form as follows: Eb = M*Ea*N, Ab = M*Aa*N, 
% Bwb = M*Bwa, Bub = M*Bua*Q, Czb = Cza*N, Dzwb = Dzwa, Dzub = Dzua*Q, 
% Cyb = P*Cya*N, Dywb = P*Dywa, Dyub = P*Dyua*Q where 'a' and 'b' refer to 
% the original and the returned realization respectively.
%
% See also HINFCD.UTIL.OREDREALSVD.  

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

    % 1. Partition the plant
    assert(nu <= size(obj.B,2), 'OREDREAL: nu larger than number of inputs.');
    assert(ny <= size(obj.C,1), 'OREDREAL: ny larger than number of outputs.');
    E = obj.E;      A = obj.A;                      Bw = obj.B(:,1:(end-nu));               Bu = obj.B(:,(end-nu+1):end);
                    Cz = obj.C(1:(end-ny),:);       Dzw = obj.D(1:(end-ny),1:(end-nu));     Dzu = obj.D(1:(end-ny),(end-nu+1):end);
                    Cy = obj.C((end-ny+1):end,:);   Dyw = obj.D((end-ny+1):end,1:(end-nu)); Dyu = obj.D((end-ny+1):end,(end-nu+1):end);
    n = size(obj.A,1);

    % 2. Determine all possible order reductions
    rr = rank(Bu); r = min([nu-rank(Dzu),rr]);
    ss = rank(Cy); s = min([ny-rank(Dyw),ss]); 

    % 3. Calculate transformation matrices
         % 3.1 P and Q -> input/output scaling
         [U,S,Q] = svd(Dzu); if size(S,1)<size(S,2); S(:,(end-r+1):end) = 0; else; S((end-r+1):end,:) = 0; end; Dzub = U*S; 
         [U,S,M] = svd(Dyw'); if size(S,1)<size(S,2); S(:,(end-s+1):end) = 0; else; S(:,(end-s+1):end) = 0; end; Dywb = (U*S)'; P = M'; 

         % 3.2 M and N -> state scaling and transformation
         Bub = Bu*Q;
         [~,~,M] = svd(Bub');
         M = [M(:,(rr+1):end) M(:,1:rr)]; BubM = Bub'*M;
         M(:,(end-rr+1):end) = M(:,(end-rr+1):end)/BubM((end-rr+1):end,(end-rr+1):end);
         M = [M(:,1:(n-rr-s)) M(:,(end-rr+1):end) M(:,(n-rr-s+1):(end-rr))]';
         
         Cyb = P*Cy; 
         [~,~,N] = svd(Cyb);
         N = [N(:,(ss+1):end) N(:,1:ss)]; CybN = Cyb*N;
         N(:,(end-ss+1):end) = N(:,(end-ss+1):end)/CybN((end-ss+1):end,(end-ss+1):end);

    % 4. Transform remaining matrices and set new realization
    Eb = M*E*N;
    Ab = M*A*N;
    Bb = M*[Bw Bub]; Bb(:,(end-rr+1):end) = 0; Bb((end-rr-s+1):(end-s),(end-rr+1):end) = eye(rr); 
    Cb = [Cz ; Cyb]*N; Cb((end-ss+1):end,1:(end-ss)) = 0; Cb((end-ss+1):end,(end-ss+1):end) = eye(ss); 
    Db = [Dzw Dzub ; Dywb P*Dyu*Q];

    Eb = roundsmall(Eb,zerotol);
    Ab = roundsmall(Ab,zerotol);
    Bb = roundsmall(Bb,zerotol);
    Cb = roundsmall(Cb,zerotol);
    Db = roundsmall(Db,zerotol);

    obj = obj.set('A',Ab,'B',Bb,'C',Cb,'D',Db,'E',Eb);

end
