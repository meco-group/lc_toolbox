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

function [s] = mergestruct(varargin)
%MERGESTRUCT merge 2 or more structures
%   Merges 2 or more structures into 1 structure where the first structure
%   has priority on the second, meaning that fields with the same name will
%   be taken from the first and not from the second. When there are more
%   than 2 structures, the 2nd has priority over the 3rd and so on.

s = struct();

for k = 1:length(varargin)
    if ~isempty(varargin{k})
        s0 = varargin{k};
        fields0 = fieldnames(s0);
        fields = fieldnames(s);

        % check which fields needed merge and which are new
        [lia,~] = ismember(fields0,fields); 

        for(j = 1:length(fields0)) % iterate over fields that are set and change them to the new value
            if(lia(j)) % merge
                if(isstruct(s.(fields0{j,1})))
                    s.(fields0{j,1}) = mergestruct(s.(fields0{j,1}),s0.(fields0{j,1})); % merge the structured fields recursively
                end
            else % add
                s.(fields0{j,1}) = s0.(fields0{j,1});
            end
        end
    end
end
end

