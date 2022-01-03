classdef printer < handle
% PRINTER Printer for generating messages in the command window
% The printer allows to easily use different styles and indentations for
% messages that are displayed to the user. 

% This file is part of hinfcd.
% Copyright (c) 2019, Laurens Jacobs, MECO Research Team @ KU Leuven. 
% 
% hinfcd is free software: you can redistribute it and/or modify it under 
% the terms of the GNU Lesser General Public License as published by the 
% Free Software Foundation, version 3.
% 
% hinfcd is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for
% more details.
% 
% You should have received a copy of the GNU Lesser General Public License
% along with hinfcd. If not, see <https://www.gnu.org/licenses/>.
    
    properties
        % indentation - The current level of indentation
        indentation = 0;
    end
    
    methods
        
        function message(obj,msg,indent)
        % MESSAGE Print a (possibly formatted) message to the command window
            if nargin<3
                indent = obj.indentation;
            end
            fprintf(obj.fitOnScreen(msg,indent));
        end
        
        function error(obj,msg,stack,indent)
        % ERROR Throws a formatted error message
        % Does not finish execution! This is done by the watchdog.
            if nargin<4
                indent = obj.indentation;
            end
            msg = obj.fitOnScreen(['Error: ' msg],indent); 
            stack = strsplit(stack,'\n');
            stack = cellfun(@(x) strrep(x,'\','\\'), stack, 'un', 0);
            stacki = strjoin(cellfun(@(x) [repmat('\t',1,indent) x], stack(1:end-1), 'un', 0),'\n');
            msg = [msg stacki '\n\n'];
            fprintf(2, msg); 
        end
        
        function warning(obj,msg,indent)
        % WARNING Throws a formatted warning message
        % Does not add the warning to the MATLAB stack! This is done by the
        % watchdog. 
            if nargin<3
                indent = obj.indentation; 
            end
            msg = ['[' 8 'Warning: ' msg ']' 8];
            obj.message(msg,indent);
        end
        
        function obj = indent(obj)
        % INDENT Increase the current indentation with 1 tab
            obj.indentation = obj.indentation+1;
        end
        
        function obj = unindent(obj)
        % UNINDENT Decrease the current indentation with 1 tab
            obj.indentation = obj.indentation-1;
            if obj.indentation<0; obj.indentation = 0; end
        end
        
        function msg = fitOnScreen(obj,msg,indent)
        % FITONSCREEN Breaks lines to fit on the screen
            if nargin<3
                indent = obj.indentation;
            end
            sizdim = matlab.desktop.commandwindow.size - 2; 
            splitlines = strsplit(sprintf(msg),{'\n'},'CollapseDelimiters',false);
            nofsp = com.mathworks.services.Prefs.getIntegerPref('CommandWindowSpacesPerTab'); 
            empty = cellfun(@isempty, splitlines); 
            fulllines = splitlines(~empty); 
            msg = cell(1,length(fulllines)); 
            for i=1:length(fulllines)
                splitwords = strsplit(fulllines{i});
                thislinelength = indent*nofsp; 
                for j=1:length(splitwords)
                    thislinelength = thislinelength + length(splitwords{j}) + 1;
                    if thislinelength <= sizdim(1)
                        if j==1
                            msg{i} = [repmat('\t',1,indent) splitwords{j}];
                        else
                            msg{i} = [msg{i} ' ' splitwords{j}];
                        end
                    else
                        thislinelength = indent*nofsp + length(splitwords{j}); 
                        msg{i} = [msg{i} '\n' repmat('\t',1,indent) splitwords{j}];
                    end
                end
            end
            splitlines(~empty) = msg; 
            msg = [strjoin(splitlines,'\n') '\n']; 
        end
        
    end
end