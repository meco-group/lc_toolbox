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

function make_example(settings,names)
%MAKE_EXAMPLE Summary of this function goes here
%   Detailed explanation goes here

% create example file
mkdir([settings.folder 'example']);
example = fopen([settings.folder 'example/example.cpp'],'w');

fprintf(example,'#include "%s.h"\n#include <iostream>\n\n',settings.filename);
fprintf(example,'int main() {\n');
fprintf(example,'\tfloat input[1];\n\n');
fprintf(example,'\t%s::%s controller;\n',settings.namespace,settings.cname);
cellfun(@(x)impulse(example,settings,x),names);
fprintf(example,'}');

fclose(example);
end

function impulse(example,settings,name)

fprintf(example,'\tstd::cout << "%s" << std::endl << "=======" << std::endl;\n',name);
fprintf(example,'\tcontroller.load(%s::%s::%s);\n',settings.namespace,settings.cname,upper(name));
fprintf(example,'\tinput[0] = 1.0f/controller.Ts();\n');
fprintf(example,'\tcontroller.update(input);\n\n');
fprintf(example,'\tinput[0] = 0.0;\n');
fprintf(example,'\tfor(int k=0;k<50;k++){\n');
fprintf(example,'\t\tstd::cout << controller.actuation()[0] << std::endl;\n');
fprintf(example,'\t\tcontroller.update(input);\n\t}\n\n');

end
