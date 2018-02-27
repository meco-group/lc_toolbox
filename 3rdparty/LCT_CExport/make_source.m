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

function make_source(settings,controllers,names,Ts)

% check controllers
nstates = max(cellfun(@(x)size(x.a,1),controllers));
nout = max(cellfun(@(x)size(x.c,1),controllers));

% create source file
source = fopen([settings.folder settings.filename '.cpp'],'w');

% includes
fprintf(source,'#include "%s.h"\n\n',settings.filename);

fprintf(source,'%s::%s::%s() :\n',settings.namespace,settings.cname,settings.cname);
fprintf(source,'\t_Ts(%d),\n',Ts);
fprintf(source,'\t_algorithm(-1)\n');
fprintf(source,'{\n\t//do nothing\n}\n\n');

% load
fprintf(source,'void %s::%s::load(int algorithm) {\n',settings.namespace,settings.cname);
fprintf(source,'\tswitch(algorithm) {\n');
cellfun(@(x)fprintf(source,'\t\tcase %s:{\n\t\t\tload_%s();\n\t\t\tbreak;}\n\n',upper(x),x),names);
fprintf(source,'\t\tdefault:{\n\t\t\treturn;}\n\t}\n\n');
fprintf(source,'\t_algorithm = algorithm;\n\treset();\n}\n\n');

% update
fprintf(source,'float* %s::%s::update(float* measurements) {\n',settings.namespace,settings.cname);
fprintf(source,'\tswitch(_algorithm) {\n');
cellfun(@(x)fprintf(source,'\t\tcase %s:{\n\t\t\tupdate_%s(measurements);\n\t\t\tbreak;}\n\n',upper(x),x),names);
fprintf(source,'\t}\n\n\treturn _y;\n}\n\n');

% actuation
fprintf(source,'float* %s::%s::actuation() {\n\treturn _y;\n}\n\n',settings.namespace,settings.cname);

% Ts
fprintf(source,'float %s::%s::Ts() {\n\treturn _Ts;\n}\n\n',settings.namespace,settings.cname);

% reset
fprintf(source,'void %s::%s::reset() {\n',settings.namespace,settings.cname);
fprintf(source,'\tfor(int k = 0; k<%d; k++)\n',nstates);
fprintf(source,'\t\t_x[k] = 0.0f;\n\n');
fprintf(source,'\tfor(int k = 0; k<%d; k++)\n',nout);
fprintf(source,'\t\t_y[k] = 0.0f;\n');
fprintf(source,'}\n\n');

% private inline functions
fprintf(source,'// private inline functions\n');
cellfun(@(x,y)print_control_law(source,settings,x,y),controllers,names);

end

function print_control_law(source,settings,controller,name)
nstates = size(controller.a,1);
nout = size(controller.c,1);

% load controller
fprintf(source,'// %s functions\n',name);
fprintf(source,'inline void %s::%s::load_%s() {\n',settings.namespace,settings.cname,name);
fprintf(source,'\t//do nothing for LTI here\n');
fprintf(source,'}\n\n');

% update controller
code_out = [vector_expand(nout,'_y'); vector_expand(nstates,'_x')];
code_x = [matrix_expand(controller.c,'x'); matrix_expand(controller.a,'x')];
code_u = [matrix_expand(controller.d,'measurements'); matrix_expand(controller.b,'measurements')];
update_code = strcat({sprintf('\t')}, code_out, {' = '}, code_x, code_u, {';'});

fprintf(source,'inline void %s::%s::update_%s(float *measurements) {\n',settings.namespace,settings.cname,name);
fprintf(source,'\tfloat x[%d];\n',nstates);
fprintf(source,'\tfor(int k = 0; k<%d; k++)\n',nstates);
fprintf(source,'\t\tx[k] = _x[k];\n\n');
cellfun(@(x)fprintf(source,'%s\n',x),update_code);
fprintf(source,'}\n\n');

end
