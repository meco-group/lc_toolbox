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

function [cf,cferr]=crest2(X,N,Fe)
% CREST2 Compute the crest factors (Kr').
% Use cf=crest2(X,N,Fe) with
%  Input argument :
%  - X  : Matrix containing the Fourier coefficients of signals
%         column by column
%  - N  : Number of time samples to be used
%  - Fe : Set of the effective harmonic numbers
%  Output argument :
%  - cf : Matrix containing two row vector with the lower and
%         upper crest factor values
% See also CREST1, SCHROEDER.
 %
% P.A.N. Guillaume - version 1 / 17 November 1990
% Copyright (c) 1990 by dept. ELEC, V.U.B.
%
% Reference :  The calculation of the error on the crest factor
% approximation is based on the inequality of Berstein.
% See :  Approximation of Functions, G.G. Lorentz, p. 39.
%
if nargin<3
   if nargin<2, x=four2ti(X); else x=four2ti(X,N); end
   cf=lnorm(x,inf)./effval(X);
else
   if isempty(N), x=four2ti(X); else x=four2ti(X,N); end
   cf=lnorm(x,inf)./effval(X,Fe);
end
if nargin<2, [N,colno]=size(x); end
[rowno,colno]=size(X);
if colno==1, dummy=X; else dummy=max(abs(X)')'; end
Ft=cumsum(ones(size(dummy)));
Ft=Ft(dummy>max(dummy)*1e-6);
clear dummy
if nargin<3, Fe=Ft; end
kmax=max(max(Ft))-1;
if N>pi*kmax
   cferr=cf/(1-pi*kmax/N);
end

