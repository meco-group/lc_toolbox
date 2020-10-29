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

classdef (InferiorClasses = {?zpk,?tf,?ss,?frd}) FRDmod < Model & frd
% Create a model based on frequency domain data.
% LCToolbox counterpart of MATLAB's \c frd().
       
    methods
        function self = FRDmod(varargin)
            if (nargin>0) && (isa(varargin{1},'AbstractLFTmod'))
                assert(nargin==3,'Conversion from LFTmod to FRDmod requires 3 arguments');
                varargin{1} = freqresp(varargin{1:3});
                varargin{4} = varargin{3}; % add frequencyunit label
                varargin{3} = 'FrequencyUnit';
            end
            
            self@Model();
            self@frd(varargin{:});
        end
        
        function nin = nin(self)
            nin = size(self,2);
        end
        
        function nout = nout(self)
            nout = size(self,1);
        end
        
        function sys = std(self)
            sys = frd(self);
        end
        
        function blk = blkdiag(varargin)
            if all(cellfun(@(x)isa(x,'FRDmod'),varargin))
                idxdisc = cellfun(@(x) x.Ts~=0, varargin);
                equalTs = cell2mat(cellfun(@(x,y) x.Ts == y.Ts, varargin(idxdisc), varargin(idxdisc), 'un', 0));
                if all(equalTs)
                    % automatic discretization of continuous time frds if appropriate
                    if isempty(equalTs)
                        Ts = 0;
                    else
                        disc = varargin(idxdisc);
                        Ts = disc{1}.Ts;
                    end
                    varargin(~idxdisc) = cellfun(@(x) setfield(x, 'Ts', Ts), varargin(~idxdisc), 'un', 0);
                end
                blk = blkdiag@frd(varargin{:});     
            else
                blk = varargin{1};
                for i = 2:nargin
                    blk = blkdiag_model(blk,varargin{i});
                end
            end
        end
        
        function varargout = subsref(self,varargin)
            varargout = {subsref@frd(self,varargin{:})};
        end
        
        function varargout = size(self,varargin)       
            varargout = {size@frd(self,varargin{:})};
        end
        
        function varargout = isempty(self,varargin)       
            varargout = {isempty@frd(self,varargin{:})};
        end
        
        function sum = plus(self,other)
            sum = plus@Model(self,other);
        end
         
        function difference = minus(self,other)
            difference = minus@Model(self,other);
        end
        
        function opposite = uminus(self)
            opposite = uminus@Model(self);
        end
         
        function product = mtimes(self,other)
            product = mtimes@Model(self,other);
        end
         
        function cat = vertcat(self,varargin)
            cat = vertcat@Model(self,varargin{:});
        end
      
        function cat = horzcat(self,varargin)
            cat = horzcat@Model(self,varargin{:});
        end

    end
    
    methods (Access=protected)
        function submod = submodel(self,idxout,idxin)
            submod = self(idxout,idxin);
        end   
    end
    
    methods (Access=public)
        function varargout = bode(varargin)
            [varargout{1:nargout}] = bode@Model(varargin{:});
        end
        
        function varargout = bodemag(varargin)
            [varargout{1:nargout}] = bodemag@Model(varargin{:});
        end
        
        function varargout = bodeplot(varargin)
            [varargout{1:nargout}] = bodeplot@Model(varargin{:});
        end
        
        function varargout = nyquist(varargin)
            [varargout{1:nargout}] = nyquist@Model(varargin{:});
        end
        
        function varargout = nyquistplot(varargin)
            [varargout{1:nargout}] = nyquistplot@Model(varargin{:});
        end
        
        function varargout = impulse(varargin)
            warning('Cannot plot impulse of FRDmod');
        end
                
        function varargout = step(varargin)
            warning('Cannot plot step of FRDmod');
        end
        
        function varargout = pzmap(varargin)
            warning('Cannot plot pzmap of FRDmod');
        end
                
        function varargout = sigma(varargin)
            [varargout{1:nargout}] = sigma@Model(varargin{:});
        end
        
        function varargout = freqresp(varargin)
            [varargout{1:nargout}] = freqresp@frd(varargin{:});
        end
        
        function varargout = evalfr(varargin)
            [varargout{1:nargout}] = evalfr@frd(varargin{:});
        end
    end
end