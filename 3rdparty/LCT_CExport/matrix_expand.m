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

function [code] = matrix_expand(X,v,varargin)
%MATRIX_EXPAND expands a matrix product to C code
%   Detailed explanation goes here

options = struct('optZ',1e-10,...
                 'optO',1e-10);
if(nargin>2)
    options = mergestruct(varargin{1},options);
end

N = size(X,1);
code = cell(N,1);
for k=1:size(X,1)
    % expand cell
    code{k,1} = row_expand(X(k,:),v,options);    
end

end

function [line] = row_expand(R,v,options)
%ROW_EXPAND expand row to corresponding C code
%   Detailed explanation goes here

line = [];
for k=1:length(R)
    % expand cell
    if((options.optZ>0)&&(abs(R(k))<options.optZ))
        % skip this entry: it is zero
    elseif((options.optO>0)&&(abs(R(k)-1)<options.optZ))
        line = [line '+' v '[' num2str(k-1) '] '];  
    elseif((options.optO>0)&&(abs(R(k)+1)<options.optZ))
        line = [line '-' v '[' num2str(k-1) '] '];         
    else
        line = [line num2str(R(k),'%+5.5e') 'f*' v '[' num2str(k-1) '] '];    
    end
end

end
