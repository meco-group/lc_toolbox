function [obj,U,V] = svdreal(obj)
% SVDREAL SVD canonical realization
% 
% b = SVDREAL(a) returns b, an equivalent realization of a, such that
% E = [I 0]
%     [0 0]
%
% [b,U,V] = OREDREALSVD(obj) returns the same as the previous 
% syntax, and furthermore returns the transformation matrices that were 
% used to obtain the canonical form as follows: Aa = U*Ab*V', Ba = U*Bb, 
% Ca = Cb*V', Ea = U*Eb*V'. 
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

    re = rank(obj.E);
    [U,S,V] = svd(obj.E);
    U = U*blkdiag(sqrt(S(1:re,1:re)),eye(size(obj.E)-re));
    V = V*blkdiag(sqrt(S(1:re,1:re)),eye(size(obj.E)-re));
    A = U\obj.A/V'; 
    B = U\obj.B;
    C = obj.C/V';
    obj = obj.set('A',A,'B',B,'C',C,'E',diag([ones(1,re) zeros(1,size(obj.E,1)-re)]));
    
end