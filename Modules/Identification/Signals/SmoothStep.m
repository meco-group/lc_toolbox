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

classdef SmoothStep < TimeSignal
    %SMOOTHSTEP A 1-dimensional smooth step signal
    %   Detailed explanation goes here
    
    properties (Access = protected)
        tmin_        % start time
        tmax_        % end time
        t0_          % starting time of the step
        riseTime_    % rise time of the smooth step
    end
    
    methods
        
        function self = SmoothStep(varargin)
            
            required = {'label','fs','t0','twindow','risetime'};
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
                if strcmp(varargin{2*i-1},'risetime'); assert(~exist('risetime','var'),'You cannot define the rise time twice.'); risetime = varargin{2*i}; end;
                if ~ismember(varargin{2*i-1},[required, optional]); warning(['I could not identify the option ''' varargin{2*i-1} ''' and I am going to ignore it.']); end;
            end
            
            self = SmoothStep1(self,label,fs,twindow,t0,risetime);
            
        end
        
        function self = SmoothStep1(self, label, fs, twindow, t0, risetime)
            
            % input parsing

            p = inputParser;
            addRequired(p,'label',@(x) assert(ischar(x),'Your label should be a string.'));
            addRequired(p,'fs',@(x) assert(isnumeric(x) && length(x)==1 && x > 0,'The sample frequency should be a positive numeric scalar.'));
            addRequired(p,'twindow',@(x) assert(isnumeric(x) && (~any(size(x) ~= [1 2]) || ~any(size(x) ~= [2 1])) && x(2) > x(1) && x(1) >= 0,'Your time window [tmin, tmax] should satisfy tmax > tmin > 0.'));
            addRequired(p,'t0',@(x) assert(isnumeric(x) && length(x)==1 && x >= 0,'t0 should be a positive numeric scalar.')); 
            addRequired(p,'risetime',@(x) assert(isnumeric(x) && length(x)==1 && x > 0,'The rise time should be a positive numeric scalar.'));
            parse(p,label,fs,twindow,t0,risetime);
            
            % assign properties
            
            self.label_ = p.Results.label;
            self.tmin_ = p.Results.twindow(1);
            self.tmax_ = p.Results.twindow(2);
            self.fs_ = p.Results.fs;
            self.t0_ = p.Results.t0;
            self.riseTime_ = p.Results.risetime;
            self.fmin_ = 0;
            self.fmax_ = fs/2;
            
            % additional 'inter-argument' checks
            
            if self.t0_ > self.tmax_ || self.tmin_ > self.t0_; error('t0 does not belong to your time window.'); end;
            if self.t0_+self.riseTime_ > self.tmax_; error('The final value is not reached in your time window.'); end;

            % generate the actual signal and its spectrum (polytraj)
            
            t = self.tmin_:1/self.fs_:self.tmax_; 
            [pos,~,~,~,px] = polytraj(1,1/self.fs_,floor(self.riseTime_*self.fs_),floor(self.riseTime_*self.fs_));
            lb = sum(t > self.t0_) - length(pos);
            lu = length(t) - length(pos) - lb;
            signal = [zeros(lu,1) ; pos ; ones(lb,1)];
            self.signal_ = struct('Time',t,'Value',signal);
            
            freq = linspace(self.fmin_,self.fmax_,100);
            PX = px .* factorial(flipud((0:length(px)-1)'));
            PX = [PX ; 0];
            val = polyval(PX,1./(freq*2*pi*1i));
            self.spectrum_ = struct('Frequency',freq','Value',val);
            
        end
        
        function disp(self)
            disp(['Smooth step signal ''' self.label_ ''' with these properties:']);
            disp(['   sample frequency:     ' num2str(self.fs_) ' Hz']);
            disp(['   time window:          [' num2str(self.tmin_) ', ' num2str(self.tmax_) '] s']);
            disp(['   start time:           ' num2str(self.t0_) ' s']);
            disp(['   rise time:            ' num2str(self.riseTime_) ' s']);
        end
        
    end
end