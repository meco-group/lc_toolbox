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

function X=time2fo(x,N)
%TIME2FO Compute the Fourier serie coefficients from
% the time domain signals.
% Use X=time2fo(x,N) with
%  Input arguments :
%  - x : Matrix containing the time signals
%  - N : Number of required time samples
%  Output argument :
%  - X : Matrix containing the Fourier coefficients of
%        the different signals column by column
% See also FOUR2TI.
 %
% P.A.N. Guillaume - version 1 / 18 November 1990
% Copyright (c) 1990 by dept. ELEC, V.U.B.
%
if nargin==1
   [rowno,colno]=size(x);
   Ndummy=rowno;
   N=2;
   while Ndummy>N, N=N*2; end
end
X=fft(x,N);
X=[X(1,:);2*X(2:floor(N/2),:)]/N;

