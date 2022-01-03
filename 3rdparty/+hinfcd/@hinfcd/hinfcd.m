classdef hinfcd
% HINFCD Tool for H-infinity controller synthesis through LMIs

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

    properties (Access=public)
        % problem - The control problem
        % The control problem contains the generalized plant, input and
        % output weighting filters and the specifications.
        % See also HINFCD.PROBLEM.
        problem
        
        % opts - Options
        % The options allow the user to specify his own preferences for
        % selecting solvers, tolerances, etc. 
        % See also HINFCD.OPTIONS.
        opts
        
        % preprocessor - The problem preprocessor
        % The problem preprocessor transforms the problem into a standard
        % formulation that allows reformulation into an SDP. 
        % See also HINFCD.PREPROCESSOR.
        preprocessor
              
        % postprocessor - The problem postprocessor
        % The problem postprocessor combines the solution of the SDP and 
        % and the information of the preprocessing phase into the
        % controller that is optimized for. 
        % See also HINFCD.POSTPROCESSOR.
        postprocessor
        
        % watchdog - Monitor of the solution progress
        % The watchdog keeps track of the evolution of the different stages 
        % of the solution, save this information in a logbook, reports to
        % the user in the command window, and is responsible for error 
        % handling.  
        % See also HINFCD.WATCHDOG.
        watchdog = hinfcd.watchdog.watchdog(); 
    end
    
    methods
        function obj = hinfcd()
        end
        
        function obj = setproblem(obj,P,Wi,Wo,specs)
        % SETPROBLEM Initializes the control problem
            obj.problem = hinfcd.problem(obj.watchdog,P,Wi,Wo,specs); 
        end
        
        obj = setoptions(obj,s);
        obj = options(obj);
        [obj,gamma] = solve(obj,options);  
    end
    
end
