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

function eff=effval(X,Fe)
%EFFVAL Calculate the effective value of all the column vectors
% present in the matrix X (Fourier coefficients).
 %
% P.A.N. Guillaume - version 1 / 19 November 1990
% Copyright (c) 1990 by dept. ELEC, V.U.B.
%
X(1,:)=sqrt(2)*X(1,:);    % Correction for the DC-component
if nargin==2, X=X(Fe(:),:); end
eff=sqrt(sum(abs(X).^2)/2);

