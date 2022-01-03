classdef preprocessor 
% PREPROCESSOR Processor formulating a problem into an SDP

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
        % prob - The hinfcd problem
        % The extended H-infinity controller synthesis problem.
        prob
        
        % opts - The options for controller synthesis
        opts
    end
    
    methods
        %% Constructor
        function obj = preprocessor(problem, options)
        % PREPROCESSOR Constructor
            
            % construct the object            
            assert(isa(problem,'hinfcd.problem'), 'Invalid problem.'); 
            obj.prob = problem; 
            
            assert(isstruct(options), 'options should be a structure.');
            obj.opts = options;
            
        end
        
        %% Select the appropriate solver
        function [K,gamma] = solve(obj)
            % SOLVE Calls the appropriate solution method
            if strcmp(obj.opts.synthesis.formulation,'descriptor')
                [K,gamma] = hinfcd.standard.projlem.projlem_d(obj.prob, obj.opts); 
            else
                [K,gamma] = hinfcd.standard.projlem.projlem_ss(obj.prob, obj.opts);
            end
        end
        
        %% Solve standard problem in descriptor form
        function [K,gamma] = solve_descriptor(obj)
            % SOLVE_DESCRIPTOR Calls the standard descriptor controller synthesis procedure 
            [K,gamma] = hinfcd.standard.projlem.projlem_d(obj.prob, obj.opts);  
        end
        
        
        %% Solve standard problem in state-space form
        function [K,gamma] = solve_statespace(obj)
            % SOLVE_STATESPACE Calls the standard state-space controller synthesis procedure
            [K,gamma] = hinfcd.standard.projlem.projlem_ss(obj.prob, obj.opts);  
        end
           
    end
    
end