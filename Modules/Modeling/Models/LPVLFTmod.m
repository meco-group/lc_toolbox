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

classdef(InferiorClasses = {?zpk,?tf,?ss,?frd}) LPVLFTmod < AbstractLFTmod & AbstractLPVmod & Model
    %LFTmod Summary of this class goes here
    %   Detailed explanation goes here
        
    methods
        function self = LPVLFTmod(M,Nu,Nl,E,parameters,varargin)
            self@AbstractLFTmod(M,Nu,Nl,E,varargin{:});
            self@AbstractLPVmod(parameters);
        end
        
        function sys = std(self)
            sys = std(gridme(self));
        end
        
        function sys = simplify(self)
            S = transpose(lft2ss(self));
            sys = LPVDSSmod(S{:},self.E,self.parameters(),self.Ts);
        end
        
        function mod = gridme(self,varargin)
            function chk = isvalidgrid(grid)
                isc = iscell(grid); 
                hastwocol = size(grid,2) == 2;
                if ~(hastwocol && isc); chk = false; return; end
                containsparam = all(cell2mat(cellfun(@ischar, grid(:,1), 'un', 0)));
                containsnumval = all(cell2mat(cellfun(@isnumeric, grid(:,2), 'un', 0)));
                chk = isc && containsparam && containsnumval;
            end
            
            if nargin>1 && isvalidgrid(varargin{1})
                args = varargin{1}(:,1)';
                grid = varargin{1}(:,2)';
            else
                [grid,args] = makegrid(self,varargin{:});
            end
            [Ms,Nus,Nls,Es] = grid_eval(self,grid,args);
            
            % Make cellgrid of LTILFTmodels
            density = cellfun(@length,grid);
            d = num2cell(density);
            if length(d) == 1, d{2} = 1; end
            C = cell(d{:});
            for k = 1:prod(density)
                C{k} = LTILFTmod(Ms(:,:,k),Nus(:,:,k),Nls(:,:,k),Es(:,:,k),self.Ts);
            end
            
            mod = Gridmod(C,[args(:),grid(:)]);
        end
        
        function mod = evalme(self,val,args)
            function p = getargument(p)
                if isa(p,'SchedulingParameter')
                    p = p.tensor_basis.arguments;
                end
            end
            
            if ~iscell(args)
                args = {args};
            end
            args = cellfun(@(x) getargument(x),args,'un',0); % behind the scenes args are always strings
            [Ms,Nus,Nls,Es] = self.eval(val,args);
            mod = LTILFTmod(squeeze(Ms),squeeze(Nus),squeeze(Nls),squeeze(Es),self.Ts);
        end
        
        function product = mtimes(self,other)
            if isnumeric(self) || isnumeric(other)
                product = mtimes@AbstractLFTmod(self,other);
            else
                product = mtimes@Model(self,other);
            end
        end
    end
end