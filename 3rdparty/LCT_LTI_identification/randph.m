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

function X=randph(X) 
%RANDPH random generation of phase between -pi and +pi. 
% Use X=randph(X) with 
%  Input argument : 
%  - X : abs(X) is a matrix containing the ampl. spectrum of the 
%        different signals column by column. 
%  Output argument : 
%  - X : Matrix containing the spectrum with random phases 
%        phase.n = phase.1 - 2p.[S{k=1,n-1} (n-k)Ak^2]. 
% 
ampl=abs(X); 
[freqno,signo]=size(ampl); 
%rand('uniform'); 
rand('seed',sum(100*clock)); 
phase = (2*rand(size(ampl))-1)*pi; 
X=ampl.*exp(j*phase); 
 
