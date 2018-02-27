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

function cost = mlfdi_res(Bn,An,freq,X,Y,sX2,sY2,cXY,cORd,fs)
%
%  cost = mlfdi_res(Bn,An,freq,X,Y,sX2,sY2,cXY,cORd,fs)
%

j=sqrt(-1);
freq=freq(:);
N=length(freq);
[x,n]=size(An);
n=n-1;

% calculation of the frequency axis
if (cORd == 'c')
   waxis = j*2*pi*freq;
elseif (cORd == 'd')
   waxis = exp(j*2*pi*freq/fs);
else
   disp('time domain is undefined; it is set to continuous time');
   cORd = 'c';
   waxis = j*2*pi*freq;
end;


% matrices P and Q are used to form the complex set of equations
P = kron(ones(1,n+1),waxis).^kron(ones(N,1),(n:-1:0));

Num = P*Bn';
Den = P*An';
SE = sqrt(sX2.*(abs(Num).^2) + sY2.*(abs(Den).^2) - 2*real(cXY.*Den.*conj(Num)) );
E = (Num.*X - Den.*Y)./SE;

cost = (norm(E)^2)/2;
