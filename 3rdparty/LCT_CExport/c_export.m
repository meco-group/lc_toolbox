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

function c_export(settings,varargin)
%C_EXPORT Summary of this function goes here
%   Detailed explanation goes here

assert(mod(length(varargin),2)==0,'expected even number of controllers and names');

% split input arguments
controllers = varargin(1:2:end);
names = varargin(2:2:end);
Ts = controllers{1}.Ts;

sprintf('Exporting %d %s to %s',length(controllers),settings.cname,settings.filename);
    
% make header
make_header(settings,controllers,names);
% make source
make_source(settings,controllers,names,Ts);

if(isfield(settings,'example') && settings.example)
  % make example
  make_example(settings,names);
  % make makefile
  make_makefile(settings);
end

end

