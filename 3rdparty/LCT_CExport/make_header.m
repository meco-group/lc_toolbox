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

function make_header(settings,controllers,names)

% check controllers
nstates = max(cellfun(@(x)size(x.a,1),controllers));
nout = max(cellfun(@(x)size(x.c,1),controllers));

% create header file
header = fopen([settings.folder settings.filename '.h'],'w');

% include guard
fprintf(header,'#ifndef %s_H\n',upper(settings.cname));
fprintf(header,'#define %s_H\n\n',upper(settings.cname));

% namespace
fprintf(header,'namespace %s {\n\n',settings.namespace);

% class
fprintf(header,'class %s\n{\n',settings.cname);

% namespace and enum
fprintf(header,'public:\n');
fprintf(header,'\ttypedef enum algorithm_t {\n');
cellfun(@(x) fprintf(header,'\t\t%s,\n',upper(x)),names);
fprintf(header,'\t} algorithm_t;\n\n');

fprintf(header,'private:\n\tconst float _Ts;\n\tint _algorithm;\n\n');
fprintf(header,'\tfloat _x[%d];\n',nstates);
fprintf(header,'\tfloat _y[%d];\n\n',nout);

% load and update per controller
cellfun(@(x) fprintf(header,'\t//%s functions\n\tinline void load_%s();\n\tinline void update_%s(float* measurements);\n\n',x,x,x),names);

fprintf(header,'public:\n');
fprintf(header,'\t%s();\n\n',settings.cname);

fprintf(header,'\tvoid load(int algorithm);\n');
fprintf(header,'\tfloat* update(float* measurements);\n\n');

fprintf(header,'\tfloat* actuation();\n');
fprintf(header,'\tfloat Ts();\n');
fprintf(header,'\tvoid reset();\n');

fprintf(header,'}; //class\n\n');
fprintf(header,'}; //namespace\n\n');
fprintf(header,'#endif //%s_H\n',upper(settings.cname));

fclose(header);
end
