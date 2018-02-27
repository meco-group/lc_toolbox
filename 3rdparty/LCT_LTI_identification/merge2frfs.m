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

function [ mergedFRF ] = merge2frfs( FDmeas1, FDmeas2 )

assert(nargin == 2, 'Incorrect number of input arguments.')
Gnonp = [FDmeas1.Response; FDmeas2.Response];
[freq, indFreq] = sort([FDmeas1.Frequency; FDmeas2.Frequency]);
[freqUniq, indFreqUniq] = unique(freq);
Gnonp = Gnonp(indFreq(indFreqUniq));
mergedFRF = struct('Frequency',freqUniq,'Response',Gnonp);

end

