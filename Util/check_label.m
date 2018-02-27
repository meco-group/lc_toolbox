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

function [match,m] = check_label(sys,io)
%CHECK_LABEL Check the consistency input or output labels of a system
switch io
    case 'i'
        list = sys.u;
    case 'o'
        list = sys.y;
end

    p = regexp(list,'\w*(?=(\w*))', 'match');

    if(size(p,1)==length(list))
        m = p{1,1}{1,1};
        match = 1; k = 2;
        while((match==1)&&(k<=length(list)))
            match = match && strcmp(m,p{k,1}{1,1});
            k = k+1;
        end
    else
        match = 0; m = 0;
    end

end

