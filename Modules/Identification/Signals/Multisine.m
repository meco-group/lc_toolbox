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

classdef  Multisine < TimeSignal
    %MULTISINE A 1-dimensional multisine signal
    %   Detailed explanation goes here
    
    properties (Access = protected)
        freqres_     % desired frequency resolution
        phases_      % phase relation of the different components (Schroeder or random)
        mstype_      % the multisine type (full, odd, odd-odd)
    end
    
    methods
        
        function self = Multisine(varargin)
            
            required = {'label','fs','fwindow'};
            optional = {'freqres','ppp','type','phases','amplitude'};
            
            calls = varargin(1:2:end);
            
            assert(mod(nargin,2) == 0,'Odd number of arguments.');
            assert(iscellstr(calls),'The names of your name/value pairs should be strings.');
            assert(~any(~ismember(required,calls)),['These arguments are required: ' strjoin(required)]);
            assert(sum(ismember({'ppp','freqres'},calls))==1,'You should provide either ''ppp'' (points per period) or ''freqres'' (the desired frequency resolution).');
            
            for i = 1:nargin/2
                if strcmp(varargin{2*i-1},'label'); assert(~exist('label','var'),'You cannot define a label twice.'); label = varargin{2*i}; end;
                if strcmp(varargin{2*i-1},'fs'); assert(~exist('fs','var'),'You cannot define a sample frequency twice.'); fs = varargin{2*i}; end;
                if strcmp(varargin{2*i-1},'fwindow'); assert(~exist('fwindow','var'),'You cannot define a frequency window twice.'); fwindow = varargin{2*i}; end;
                if strcmp(varargin{2*i-1},'freqres'); assert(~exist('freqres','var'),'You cannot define a frequency resolution twice.'); freqres = varargin{2*i}; end;
                if strcmp(varargin{2*i-1},'ppp'); assert(~exist('ppp','var'),'You cannot define the points per period twice.'); ppp = varargin{2*i}; end;
                if strcmp(varargin{2*i-1},'type'); assert(~exist('type','var'),'You cannot define a multisine type twice.'); type = varargin{2*i}; end;
                if strcmp(varargin{2*i-1},'phases'); assert(~exist('phases','var'),'You cannot define the phases twice.'); phases = varargin{2*i}; end;
                if strcmp(varargin{2*i-1},'amplitude'); assert(~exist('amplitude','var'),'You cannot define an amplitude spectrum twice.'); amplitude = varargin{2*i}; end;
                if ~ismember(varargin{2*i-1},[required, optional]); warning(['I could not identify the option ''' varargin{2*i-1} ''' and I am going to ignore it.']); end;
            end
            
            if exist('ppp','var')
                freqres = fs/ppp;
            end
            
            if ~exist('phases','var'); phases = 'random'; end;
            if ~exist('type','var'); type = 'full'; end;
            if ~exist('amplitude','var'); amplitude = tf(1,1); end;
            
            self = Multisine1(self,label,fs,fwindow,freqres,type,phases,amplitude);
            
        end
        
        function freq = excf(self)
            freq = self.spectrum_.Excited;
        end
        
        function out = freqres(self)
            out = self.freqres_;
        end
        
        function disp(self)
            disp(['Multisine ''' self.label_ ''' with these properties:']);
            disp(['   sample frequency:     ' num2str(self.fs_) ' Hz']);
            disp(['   frequency range:      [' num2str(self.fmin_) ', ' num2str(self.fmax_) '] Hz']);
            disp(['   frequency resolution: ' num2str(self.freqres_) ' Hz']);
            
            if strcmp(self.phases_,'schroeder'); ph = 'Schroeder'; else ph = 'random'; end; 
            if strcmp(self.mstype_,'oddodd'); mst = 'odd-odd'; else mst = self.mstype_; end;
            
            disp(['   ' mst ' multisine with ' ph ' phases']);
        end
        
    end
    
    methods (Access = private)
         function self = Multisine1(self,label, fs, fwindow, freqres, type, phases, amplitude)
            
            % input parsing

            p = inputParser;
            addRequired(p,'label',@(x) assert(ischar(x),'Your label should be a string.'));
            addRequired(p,'fs',@(x) assert(isnumeric(x) && length(x)==1 && x > 0,'The sample frequency should be a positive numeric scalar.'));
            addRequired(p,'fwindow',@(x) assert(isnumeric(x) && (~any(size(x) ~= [1 2]) || ~any(size(x) ~= [2 1])) && x(2) > x(1) && x(1) >= 0,'Your frequency window [fmin, fmax] should satisfy fmax > fmin >= 0.'));
            addRequired(p,'freqres',@(x) assert(isnumeric(x) && length(x)==1 && x > 0,'The frequency resolution should be a positive numeric scalar.'));
            addRequired(p,'type',@(x) assert(any(validatestring(x,{'full','odd','oddodd'})),'The only valid options for ''type'' are ''full'', ''odd'' or ''oddodd''.'));
            addRequired(p,'phases',@(x) assert(any(validatestring(x,{'random','schroeder'})),'The only valid options for ''phases'' are ''random'' and ''schroeder''.'));
            addRequired(p,'amplitude',@(x) assert(isa(x,'numlti'),'The amplitude spectrum should be defined by an LTI system (numlti).'));
            parse(p,label,fs,fwindow,freqres,type,phases,amplitude);
            
            % assign properties
            
            self.label_ = p.Results.label;
            self.fmin_ = p.Results.fwindow(1);
            self.fmax_ = p.Results.fwindow(2);
            self.fs_ = p.Results.fs;
            self.phases_ = p.Results.phases;
            self.mstype_ = p.Results.type;
            
            nrofsamp = ceil(self.fs_/freqres);
            self.freqres_ = self.fs_/nrofsamp;
            self.fmin_ = ceil(self.fmin_/self.freqres_)*self.freqres_;
            self.fmax_ = round(self.fmax_/self.freqres_)*self.freqres_;
            
            % additional 'inter-argument' checks
            
            if self.fmax_ > self.fs_/2-self.freqres_; error('The excited frequency range cannot exceed the Nyquistfrequency.'); end;

            % generate the actual signal and its spectrum (musin_type)
            
            if strcmp(self.phases_,'schroeder'); initial = 's';
            else if(strcmp(self.phases_,'random')); initial = 'r';
                 else error('I ended up in an undefined state. Shouldn''t happen.');
                 end
            end 
            
            if strcmp(self.mstype_,'oddodd'); tp = 'O';
            else if strcmp(self.mstype_,'odd'); tp = 'o';
                 else if strcmp(self.mstype_,'full'); tp = 'f';
                      else error('I ended up in an undefined state. Shouldn''t happen.');
                      end
                 end
            end
            
            cp = 'c'; 
           
            [Bn,An] = tfdata(p.Results.amplitude); Bn = Bn{:}; An = An{:};
            
            [x,X,f,~,ft] = musin_type(self.fs_,self.fmin_,self.fmax_,nrofsamp,initial,cp,tp,Bn,An);
            x = x/max(abs(x)); 
            
            switch self.mstype_
                case 'oddodd'
                    nminexc = floor(ceil(self.fmin_/self.freqres_)/4)*4+2;
                    nmaxexc = round(self.fmax_/self.freqres_);
                    fexc = ft(nminexc:4:nmaxexc);
                    f = ft(nminexc:nmaxexc);
                case 'odd'
                    nminexc = floor(ceil(self.fmin_/self.freqres_)/2)*2+2;
                    nmaxexc = round(self.fmax_/self.freqres_);
                    fexc = ft(nminexc:2:nmaxexc);
                    f = ft(nminexc:nmaxexc);
                case 'full'
                    fexc = f;
            end
            
            t = 0:1/self.fs_:(nrofsamp-1)/self.fs_;
            self.signal_ = struct('Time',t','Value',x);
            self.spectrum_ = struct('Frequency',f,'Value',X,'Excited',fexc);
            
         end
        
    end
end

