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

classdef Parabola < MonomialSignal
    %RAMP A 1-dimensional step signal
    %   Detailed explanation goes here
    
    methods
        function self = Parabola(varargin)
            
            required = {'label','fs','t0','twindow'};
            optional = {};
            
            calls = varargin(1:2:end);
            
            assert(mod(nargin,2) == 0,'Odd number of arguments.');
            assert(iscellstr(calls),'The names of your name/value pairs should be strings.');
            assert(~any(~ismember(required,calls)),['These arguments are required: ' strjoin(required)]);
            
            for i = 1:nargin/2
                if strcmp(varargin{2*i-1},'label'); assert(~exist('label','var'),'You cannot define a label twice.'); label = varargin{2*i}; end;
                if strcmp(varargin{2*i-1},'twindow'); assert(~exist('twindow','var'),'You cannot define a time window twice.'); twindow = varargin{2*i}; end;
                if strcmp(varargin{2*i-1},'t0'); assert(~exist('t0','var'),'You cannot t0 twice.'); t0 = varargin{2*i}; end;
                if strcmp(varargin{2*i-1},'fs'); assert(~exist('fs','var'),'You cannot define the sample frequency twice.'); fs = varargin{2*i}; end;
                if ~ismember(varargin{2*i-1},[required, optional]); warning(['I could not identify the option ''' varargin{2*i-1} ''' and I am going to ignore it.']); end;
            end
            
            self = self@MonomialSignal('degree',2,'label',label,'twindow',twindow,'t0',t0,'fs',fs);
            
        end
        
        function disp(self)
            disp(['Parabolic signal ''' self.label ''' with these properties:']);
            disp(['   sample frequency:     ' num2str(self.fs_) ' Hz']);
            disp(['   time window:          [' num2str(self.tmin_) ', ' num2str(self.tmax_) '] s']);
            disp(['   start time:           ' num2str(self.t0_) ' s']);
        end
        
    end
end