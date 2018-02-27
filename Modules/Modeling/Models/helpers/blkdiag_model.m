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

function blk = blkdiag_model(self,other)
%BLKDIAG_LFT_FRD Blkdiag LFTmod and FRDmod
%   Compute block-diagonal of an Umod, LFTmod and FRDmod. The output is of the
%   type FRDmod, in case of Umod and FRDmod, it returns with the nominal model alongwith uncertain channel.
%

if isempty(self)
    blk = other;
elseif isempty(other)
    blk = self;
else
    if isnumeric(self), self = SSmod(self); end
    if isnumeric(other), other = SSmod(other); end
    if isa(self,'AbstractLFTmod')
        if isa(other,'AbstractLFTmod')
            blk = blkdiag(self,other);
        elseif isa(other,'FRDmod')
            self = FRDmod(self,other.Frequency,other.FrequencyUnit);
            self.Ts = other.Ts;
            blk = blkdiag(self,other);
        elseif isa(other,'Gridmod')
            blk = blkdiag(self,other);
        elseif isa(other,'Umod')
            blk = blkdiag(Umod(self),other);
        elseif isa(other,'ODEmod')
            blk = blkdiag(ODEmod(self),other);
        else
            error('Arguments must be of type LFTmod, FRDmod, Umod or ODEmod');
        end
    elseif isa(self,'FRDmod')
        if isa(other,'AbstractLFTmod')
            other = FRDmod(other,self.Frequency,self.FrequencyUnit);
            other.Ts = self.Ts;
            blk = blkdiag(self,other);
        elseif isa(other,'FRDmod')
            blk = blkdiag(self,other);
        elseif isa(other,'Umod')
            blk = blkdiag(Umod(self),other);
        elseif isa(other,'Gridmod')
            blk = blkdiag(self,other);            
        else
            error('Arguments must be of type LFTmod, FRDmod or Umod');
        end
    elseif isa(self,'Umod')
        if isa(other,'AbstractLFTmod')
            blk = blkdiag(self,Umod(other));
        elseif isa(other,'FRDmod')
            blk = blkdiag(self,Umod(other));
        elseif isa(other,'Umod')
            blk = blkdiag(self,other);
        elseif isa(other,'Gridmod')
            blk = blkdiag(self,Umod(other));
        elseif isa(other,'ODEmod')
            blk = blkdiag(self,Umod(other));
        else
            error('Arguments must be of type LFTmod, FRDmod, Umod or ODEmod');
        end
    elseif isa(self,'Gridmod')
        if isa(other,'Umod')
            blk = blkdiag(Umod(self),other);
        else
            blk = blkdiag(self,other);
        end
    elseif isa(self,'ODEmod')
        if isa(other,'AbstractLFTmod')
            blk = blkdiag(self,ODEmod(other));
        elseif isa(other,'Umod')
            blk = blkdiag(Umod(self),other);
        elseif isa(other,'ODEmod')
            blk = blkdiag(self,other);
        elseif isa(other,'Gridmod')
            blk = blkdiag(self,other);
        else
            error('Arguments must be of type LFTmod, Umod or ODEmod');
        end
    end
    
    blk.name = [self.name '-' other.name];
end
end