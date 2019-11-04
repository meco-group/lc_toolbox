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

function [args,labels] = stdargs(args,filter)
%STDARGS Convert arguments to standard form
%   Converts all models in the list to the standard form
%   Inputs
%   * args: argument list to be converted
%   * filter: remove arguments that do not fit the filter

if nargin == 1, filter = 'none'; end

% convert all systems to models
m = find(cellfun(@(x)isa(x,'AbstractSystem'),args));
models = cellfun(@(x)content(model(x)),args(m),'un',0);
for j = length(m):-1:1
    i = m(j);
    if (i+1) <= length(args) && ischar(args{i+1})
        models{j} = reshape([models{j};repmat(args(i+1),size(models{j}))],1,[]);
        args = [args(1:i-1),models{j},args(i+2:end)];
    else
        models{j} = reshape(models{j},1,[]);
        args = [args(1:i-1),models{j},args(i+1:end)];
    end
end

% filter models
switch(filter)
    case 'time'
        m = find(cellfun(@(x)isa(x,'FRDmod'),args));
        for i = fliplr(m)
            args(i) = [];
            if length(args) >= i && ischar(args{i}), args(i) = []; end
        end
    case 'freq'
        m = find(cellfun(@(x)isa(x,'ODEmod'),args));
        for i = fliplr(m)
            args(i) = [];
            if length(args) >= i && ischar(args{i}), args(i) = []; end
        end
    case 'linmod'
        m = find(cellfun(@(x)~isa(x,'AbstractLFTmod'),args));
        for i = fliplr(m)
            args(i) = [];
            if length(args) >= i && ischar(args{i}), args(i) = []; end
        end
end

% convert all models
m = cellfun(@(x)isa(x,'Model'),args);
labels = cellfun(@(x)x.name,args(m),'un',0);
args(m) = cellfun(@std,args(m),'un',0);

end

