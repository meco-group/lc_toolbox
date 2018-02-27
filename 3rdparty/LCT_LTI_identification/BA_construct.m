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

function [Bn,An]=BA_construct(y,nh,nl,M_mh,M_ml)

M_ml = M_ml(:);
M_mh = M_mh(:);

nino=length(M_ml(:));

maxB = max([M_mh;nh]);

% storing the solution in matrices An and Bn
An = [1 y(1:nh-nl)' zeros(1,nl)];
Bn = zeros(nino,maxB+1);
index_count = nh-nl+1;
for (i=1:nino)
   yy = y(index_count:index_count + M_mh(i)-M_ml(i))';
%   Bn(i,maxB-M_mh(i):maxB-M_ml(i))=yy;
   Bn(i,maxB+1-M_mh(i):maxB+1-M_ml(i))=yy;
 
   index_count = index_count + M_mh(i)-M_ml(i) + 1;
end;
