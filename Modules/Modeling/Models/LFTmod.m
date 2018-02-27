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

function mod = LFTmod(M,Nu,Nl,E,varargin)
%LFTMOD Unified LFTmod construction
%   This function determines the simplest model structure for the given
%   input. mod will either be of type LPVLFTmod or LTILFTmod.

if ~all(size(M)==[4,4])
    n = size(M) - (size(Nu) + size(E) + size(Nl));
    dims = [size(Nu)',size(E)',n',size(Nl)'];
    M = mat2cell(M,dims(1,:),dims(2,:));
end
if iscell(M), checkM = all(cellfun(@isnumeric,M(:)));
else checkM = isnumeric(M); end
if checkM && isnumeric(Nu) && isnumeric(Nl) && isnumeric(E)
    mod = LTILFTmod(M,Nu,Nl,E,varargin{:});
else
    mod = LPVLFTmod(M,Nu,Nl,E,varargin{:});
end

end