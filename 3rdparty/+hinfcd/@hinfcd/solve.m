function [K,gamma] = solve(obj,options) 
% SOLVE Solves the control problem
    
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

    import hinfcd.util.*;

    if nargin>1 && isstruct(options)
        obj = obj.setoptions(options); 
    else
        obj = obj.setoptions(); 
    end
    
    % 0. Start
    obj.watchdog.info(bold('THIS IS HINFCD, VERSION 1.0 (May 2020)'));
    
    % 1. Formulate the problem as an SDP
    obj.watchdog.info('STARTED PREPROCESSING');
    obj.watchdog.startSubTask();
        obj.watchdog.info('Formulating the SDP...');
        obj.preprocessor = hinfcd.preprocessor(obj.problem, obj.opts);
    obj.watchdog.endSubTask();
    
    % 2. Solve the SDP
    obj.watchdog.info('STARTED SOLVING');
    obj.watchdog.startSubTask();
        obj.watchdog.info('Solving the SDP...');
        [K,gamma]= obj.preprocessor.solve();
    obj.watchdog.endSubTask();
    
    % 3. Postprocess the solution
    obj.watchdog.info('STARTED POSTPROCESSING...');
    obj.postprocessor = hinfcd.postprocessor(obj.problem, K, gamma);
    obj.watchdog.startSubTask();
        obj.watchdog.info('Recombining the controller with the unstable and improper modes of the weighting filters...'); 
        obj.postprocessor = obj.postprocessor.recombine();
        obj.postprocessor = obj.postprocessor.elimstaticmodes(); 
        obj.postprocessor = obj.postprocessor.compensateDyu();
    obj.watchdog.endSubTask(); 
    
    % 4. Return
    K = obj.postprocessor.K; 
    
end