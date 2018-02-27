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

function [ output_args ] = make_makefile(settings)
%MAKE_MAKEFILE Summary of this function goes here
%   Detailed explanation goes here

% create makefile
makefile = fopen([settings.folder 'example/makefile'],'w');

fprintf(makefile,'example: example.cpp ../%s.cpp\n',settings.filename);
fprintf(makefile,'\tg++ -o example example.cpp ../%s.cpp -I ..\n',settings.filename);

fclose(makefile);
end

