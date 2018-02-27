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

function X=schroed(X)
%SCHROED Schroeder phase coding to obtain low crest factors.
% Use X=schroed(X) with
%  Input argument :
%  - X : abs(X) is a matrix containing the ampl. spectrum of the
%        different signals column by column.
%  Output argument :
%  - X : Matrix containing the spectrum with shroeder phases
%        phase.n = phase.1 - 2p.[S{k=1,n-1} (n-k)Ak^2].
% See also CREST1, CREST2, CFMINIMAX.
 %
% P.A.N. Guillaume - version 1 / 19 May 1990
% Copyright (c) 1990 by dept. ELEC, V.U.B.
%
ampl=abs(X);
[freqno,signo]=size(ampl);
phase=zeros(size(ampl));
amplnorm=ampl./(ones(size(freqno,1))*sqrt(sum(ampl.^2)));
amplnorm=2*pi*amplnorm.^2;
for i=3:freqno,
   phase(i,:)=phase(i-1,:)-sum(amplnorm(1:i-1,:));
end
X=ampl.*exp(j*phase);

