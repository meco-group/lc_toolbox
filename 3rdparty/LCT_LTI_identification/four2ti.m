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

function x=four2ti(X,N)
%FOUR2TI Compute the time domain signals from the Fourier
% serie coefficients.
% Use x=four2ti(X,N) with
%  Input arguments :
%  - X : Matrix containing the Fourier coefficients of
%        the different signals column by column
%  - N : Number of required time samples
%  Output argument :
%  - x : Matrix containing the time signals
% See also TIME2FO.
 %
% P.A.N. Guillaume - version 1 / 17 November 1990
% Copyright (c) 1990 by dept. ELEC, V.U.B.
%
if nargin==1
   [rowno,colno]=size(X);
   Ndummy=rowno*2;
   N=2;
   while Ndummy>N, N=N*2; end
end
x=N*real(ifft(X,N));

