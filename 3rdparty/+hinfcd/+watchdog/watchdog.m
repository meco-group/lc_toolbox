classdef watchdog < handle
% WATCHDOG Monitor for the solution process of the control problem
% The watchdog keeps track of the evolution of the different stages of the
% solution procedure, save this information in a logbook, reports to
% the user in the command window, and is responsible for error 
% handling.  

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
        % logbook - The logbook of the solution procedure
        % Every event that is reported by the different phases during the
        % solution procedure are saved in this logbook. 
        logbook
        
        % printer - Printer generating messages to the user
        % The output to the terminal/command window is handled by a
        % separate printer to facilitate its lay-out. 
        printer
    end
    
    methods
        
        function obj = watchdog()
        % WATCHDOG Constructor
            obj.printer = hinfcd.watchdog.printer(); 
        end
        
        function writeToLog(obj,type,msg)
        % WRITETOLOG Adds a new entry to the logbook
            import hinfcd.util.*;
            if ~isempty(msg)
                obj.logbook(end+1).time = datestr(now);
                obj.logbook(end).tag = type;
                obj.logbook(end).msg = removehtml(msg);
            end
        end
        
        function error(obj,msg,internal)
        % ERROR Throws a formatted error
            if nargin<3
                internal = 0;
            end
            stack = evalin('caller','evalc(''dbstack'')');
            if internal
                split = strsplit(stack,'\n');
                stack = strjoin(split(1+internal:end),'\n');
                stack = ['>' stack(2:end)];
            end
            obj.printer.error(msg,stack);
            obj.writeToLog('ERROR',[msg '\n' stack]);
            err = struct(); 
            err.message = 'hinfcd detected an error in the execution of your command. See the former error message or the logbook for more information.'; 
            err.identifier = 'hinfcd:watchdogError'; 
            err.stack.file = '';
            err.stack.name = 'hinfcd';
            err.stack.line = 1;
            error(err); 
        end
        
        function warning(obj,msg,internal)
        % WARNING Throws a formatted warning
            if nargin<3
                internal = 0;
            end
            stack = evalin('caller','evalc(''dbstack'')');
            if internal
                split = strsplit(stack,'\n');
                stack = strjoin(split(1+internal:end),'\n');
                stack = ['>' stack(2:end)];
            end
            obj.printer.warning(msg);
            obj.writeToLog('WARNING',[msg '\n' stack]);
            w = warning(); 
            warning('off','all');
            warning(msg);
            warning(w); 
        end
        
        function info(obj,msg)
        % INFO Displays a formatted message
            obj.writeToLog('INFO',msg);
            obj.printer.message(msg);
        end

        function assertCond(obj,condition,msg)
        % ASSERTCOND Asserts a logical condition 
            if nargin<3
                msg = 'Assertion failed.';
            end
            if ~condition
                obj.error(msg,1);
            end
        end
        
        function assertType(obj,object,type,msg)
        % ASSERTTYPE Asserts a specific object type
            if nargin<4
                msg = ['Object is not of type ''' type '''.'];
            end
            if ~isa(object,type)
                obj.error(msg,1);
            end
        end
        
        function assertField(obj,object,field,msg)
        % ASSERTFIELD Asserts the existence of a certain field in a structure (array)
            if nargin<4
                msg = ['A field with name ''' field ''' is missing.'];
            end
            if ~isstruct(object) || ~isfield(object,field)
                obj.error(msg,1);
            end
        end
                
        function assertNonEmpty(obj,object,msg)
        % ASSERTNONEMPTY Asserts that an object is not empty
            if nargin<3
                msg = 'An object or value was not expected to be empty.';
            end
            if isempty(object)
                obj.error(msg,1); 
            end
        end
        
        function warnCond(obj,condition,msg)
        % WARNCOND Throws a warning 
            if ~condition
                obj.warning(msg,1);
            end
        end
        
        function startSubTask(obj)
        % STARTSUBTASK Indicates that the execution of a subtask is initialized
            obj.printer.indent();
        end
        
        function endSubTask(obj)
        % ENDSUBTASK Indicates that a subtask is finished
            obj.printer.unindent();
        end
        
        function varargout = executeSilently(obj, cmd, msg)
        % EXECUTESILENTLY Executes MATLAB functions while suppressing warnings and information
            try
                [log,varargout{1:nargout}] = evalc(['evalin(''caller'',''' cmd ''')']);
                obj.writeToLog('INFO', log);
                if ~isempty(strfind(log,'Warning:'))
                    obj.writeToLog('WARNING','Silenced execution threw (a) warning(s). Check the former INFO message.'); 
                end
            catch hiddenerr
                obj.writeToLog('INFO', getReport(hiddenerr));
                obj.error(msg,2);
            end
        end
            
    end
end