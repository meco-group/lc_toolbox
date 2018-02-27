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

function Lp=lnorm(A,p)
%LNORM Calculate the vector Lp-norm of all the column vectors
% present in the matrix A.
% See also NORM.
 %
% P.A.N. Guillaume - version 1 / 18 November 1990
% Copyright (c) 1990 by dept. ELEC, V.U.B.
%
if nargin==1, p=2; end
[rowno,colno]=size(A);
for i=1:colno, Lp(i)=norm(A(:,i),p); end
if p~=inf, Lp=Lp/(rowno^(1/p)); end

