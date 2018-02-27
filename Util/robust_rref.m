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

function [RR] = robust_rref(M)
%ROBUST_RREF Summary of this function goes here
%   Detailed explanation goes here

cm = cond(M);
if(cm>1e10)
    warning(['M is badly conditioned: ' num2str(cm)])
end

[R,ja] = rref(M,1e-16);
r = length(ja);

% tol = max(size(M)) * eps(norm(M));

jb = 1:size(M,2);
jb(ja) = [];

rankcheck = rank(M,1e-16) - r;
if(rankcheck ~= 0)
    warning(['rank error: rank(M) = ' num2str(rank(M)) ', length(r) = ' num2str(r)]);
end
    
rankcheck = rank(M(1:r,ja)) - r;
% tol = r * eps(norm(M(1:r,ja)));

if(rankcheck ~= 0)
    warning(['rank error: rankcheck = ' num2str(rankcheck)]);
end

S = M(1:r,ja)\M(1:r,jb);
% S = linsolve(M(1:r,ja),M(1:r,jb));

RR = zeros(size(M));
RR(1:r,ja) = eye(r,r);
RR(1:r,jb) = S;

d = R - RR;
if(max(abs(d(:)))>1e-14)
    warning(['maximum tolerance surpassed' num2str(max(abs(d(:))))])
end

end

