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

function sys = fromstd(sys,varargin)
% Creates an LCToolbox model from a standard MATLAB \c numlti model.
%
% Parameters:
%  sys : the standard MATLAB model @type numlti
%  varargin : may contain a char with the name of the model
% 
% Return values:
%  sys : the equivalent LCToolbox model @type Model

c = cellfun(@ischar,varargin);
assert(sum(c)<=1,'Only one string argument allowed');
if any(c)
    name = varargin{c};
    varargin(c) = [];
end
if ~isa(sys,'Model')
    if ndims(sys) <= 2
        if isa(sys,'frd')
            p = properties(sys);
            sys_ = FRDmod(sys.ResponseData,sys.Frequency);
            for i = 1:length(p)
                sys_.(p{i}) = sys.(p{i});
            end
            sys = sys_;
        elseif isnumeric(sys)
            sys = SSmod(sys,varargin{:});
        else
            [A,B,C,D,E,Ts] = dssdata(sys);
            sys = DSSmod(A,B,C,D,E,Ts,varargin{:});
        end
    else
         argout = cell(ndims(sys),1);
         argout = size(sys);
         siz = argout(3:end);

         idcs = cell(1,ndims(sys)-2);
         grid = cell(siz);
         for i = 1:prod(siz)
            [idcs{:}] = ind2sub(siz,i);
            eval(['grid{idcs{:}} = fromstd(sys(:,:,' strjoin(cellfun(@(x) num2str(x),idcs,'un',0),',') '));']);
         end
         
         params(:,1) = fieldnames(sys.SamplingGrid);
         params(:,2) = struct2cell(sys.SamplingGrid);
         for i = 1:length(params(:,1))
             idcs(:) = {1}; idcs{i,1} = ':';
             params{i,2} = subsref(params{i,2},substruct('()',idcs));
             params{i,2} = reshape(params{i,2},1,length(params{i,2}));
         end
         sys = Gridmod(grid,params);
    end
end
if any(c)
    if ~isempty(sys.name) && ~strcmp(sys.name,name)
        warning(['Overwriting model name: ',sys.name, ' by ' ,name]);
    end
    sys.name = name;
end
end

