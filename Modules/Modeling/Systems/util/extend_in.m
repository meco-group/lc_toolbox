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

function varargout = extend_in(in,varargin)
    if nargin > 2
        varargout = cellfun(@(x)extend_in(in,x),varargin,'UniformOutput',false);
    else
        self = varargin{1};
        M_ = self.M;            

        [loca,locin] = ismember(self.in,in);
        B_ = zeros(size(M_,1),length(in));
        try
            B_(:,locin) = M_(:,self.nupper+self.nx+(1:self.nu));
        catch
            keyboard
        end
        M_ = [M_(:,1:(self.nupper+self.nx)), B_,...
            M_(:,self.nupper+self.nx+self.nu+(1:self.nlower))];

        varargout{1} = LFTsys(M_,self.Nu,self.Nl,self.E,self.Ts,in,self.out);
    end
end
