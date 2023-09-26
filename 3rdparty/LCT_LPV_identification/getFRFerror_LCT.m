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

function [ errorFRF ] = getFRFerror_LCT(freqLines, OmegaConc, FRFmConc, pGridIdent, schParam, model )

FRFerror = zeros(numel(FRFmConc),1);
for k = 1:numel(pGridIdent)
    mdlEval = model.evalme(pGridIdent(k),{schParam}); 
    FRFerror(((freqLines(k)-1)*nout(mdlEval)*nin(mdlEval)+1):((freqLines(k+1)-1)*nout(mdlEval)*nin(mdlEval)),1) = vec(freqresp(mdlEval,OmegaConc((freqLines(k):(freqLines(k+1)-1)),1))) - FRFmConc(((freqLines(k)-1)*nout(mdlEval)*nin(mdlEval)+1):((freqLines(k+1)-1)*bout(mdlEval)*nin(mdlEval)),1); 
end
errorFRF = [real(FRFerror); imag(FRFerror)];

end

