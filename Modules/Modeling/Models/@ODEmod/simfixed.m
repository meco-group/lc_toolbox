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

function varargout = simfixed(type,varargin)
    function [ys,ts] = singlesim(type,self,index,varargin)
        uv = zeros(self.nin,1); uv(index) = 1;
        switch(type)
            case 'step'
                u = @(t) uv;
            case 'impulse'
                u = @(t) uv*(100)*(t<=(1/100));
        end                
        [ys,ts,~] = sim(self,u,varargin{:},zeros(self.nx,1));
    end
    
    if nargout == 0
        if isnumeric(varargin{end})
            t0 = varargin{end};
            varargin(end) = [];
        end
        % Separate arguments
        isode = cellfun(@(x)isa(x,'ODEmod'),varargin);
        islti = cellfun(@(x)isa(x,'numlti'),varargin);
        
        ismod = find(isode | islti);
        ismode = [ismod(2:end)-1,length(varargin)];
        args = arrayfun(@(i1,i2) varargin(i1:i2),ismod,ismode,'un',0);

        if ~exist('t0','var')
            if any(islti)
                args0 = args{find(islti,1)};
                switch(type)
                    case 'step'
                        [~,t0] = step(args0{:});
                    case 'impulse'
                        t0 = []; %[~,t0] = impulse(args0{:});
                end
            else
                t0 = [];
            end
        end

        switch(type)
            case 'step'
                [y,t] = cellfun(@(x) step(x{:},t0),args,'un',0);
            case 'impulse'
                [y,t] = cellfun(@(x) impulse(x{:},t0),args,'un',0);
        end

        for i = 1:size(varargin{1},1)
            for j = 1:size(varargin{1},2)
                subplot(size(varargin{1},1),size(varargin{1},2),(i-1)*size(varargin{1},2)+j);
                nargs = cellfun(@(w,z,tau) [{tau,z(:,i,j)},w(2:end)],args,y,t,'un',0);
                hargs = horzcat(nargs{:});
                plot(hargs{:});
            end
        end
    else
        if strcmp(type,'impulse')
            warning('The impulse response of an ODE model is only a (rough) approximation of the actual response to a dirac impulse.');
        end

        self = varargin{1};
        varargin(1) = [];
        [y,t] = singlesim(type,self,1,varargin{:});
        isnum = cellfun(@isnumeric,varargin);
        varargin{isnum} = t;

        y2 = arrayfun(@(i) singlesim(type,self,i,varargin{:}),2:self.nin,'un',0);
        y = cat(3,y,y2{:});

        switch nargout
            case 1
                varargout = {y};
            case 2
                varargout = {y,t};
        end
    end
end


