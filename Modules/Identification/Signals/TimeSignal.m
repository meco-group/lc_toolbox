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

classdef TimeSignal
    %TIMESIGNAL A 1-dimensional time signal
    %   Detailed explanation goes here
    
    properties %(Access = protected)
        fmin_        % minimal frequency
        fmax_        % maximal frequency
        fs_          % sample frequency
        signal_      % the signal itself
        spectrum_    % the signal's continuous Fourier transform
        label_       % assign a label to the signal
        periodicity_ % whether the signal is periodic or not
    end
    
    methods
        function out = fs(self) 
            out = self.fs_; 
        end
        function out = label(self)
            out = self.label_; 
        end
        function out = type(self)
            out = class(self); 
        end
        function out = fmin(self) 
            out = self.fmin_;
        end
        function out = fmax(self) 
            out = self.fmax_; 
        end
        function out = length(self)
            out = length(self.signal); 
        end
        
        function [val,time] = signal(self,index)
            if nargin > 1
                time = self.signal_.Time(index);
                val = self.signal_.Value(index);
            else
                time = self.signal_.Time;
                val = self.signal_.Value;
            end
        end
        
        function [val,freq] = spectrum(self,varargin)
            
            p = inputParser;
            addOptional(p,'ft','dft',@(x) any(validatestring(x,{'cft','dft'})));
            parse(p,varargin{:});
            
            switch p.Results.ft
                case 'cft'
                    if isa(self,'Multisine'); warning('CFT not (yet) available for multisines. Try the ''dft'' option.'); freq = []; val = []; return; end;
                    freq = self.spectrum_.Frequency;
                    val = self.spectrum_.Value;
                case 'dft' 
                    df = self.fs_/length(self.signal_.Time);
                    nmax = floor(self.fmax_/df);
                    nmin = ceil(self.fmin_/df);
                    freq = df * (nmin:nmax)';
                    val = fft(self.signal_.Value);
                    val = val(nmin+1:nmax+1)/length(val)*2;
                otherwise
                    error('I ended up in an undefined state. Shouldn''t happen.');
            end
        end
        
        function plotSignal(self)
                figure; plot(self.signal_.Time, self.signal_.Value);
                title(['Time signal: ''' self.label_ '''']);
                xlabel('time (s)');
                ylabel('signal');
        end
            
        function plotSpectrum(self,varargin)
            
            p = inputParser;
            addOptional(p,'ft','cft',@(x) any(validatestring(x,{'cft','dft'})));
            parse(p,varargin{:});
            
            switch p.Results.ft
                case 'cft'
                    if isa(self,'Multisine'); warning('CFT not (yet) available for multisines. Try the ''dft'' option.'); return; end;
                    [val,freq] = spectrum(self,'cft');
					figure;
                    subplot(211);
                    semilogx(freq, db(val));
                    ylabel('magnitude [dB]'), title(['Spectrum (CFT) of time signal ''' self.label '''']);
					axis tight;
                    subplot(212);
                    semilogx(freq, 180/pi*unwrap(angle(val)));
                    ylabel('phase [�]'), xlabel('f  [Hz]');
                    axis tight;
                case 'dft' 
                    [val,freq] = spectrum(self,'dft');
					figure;
                    subplot(211);
                    stem(freq, db(val));
                    set(gca,'xscal','log');
                    ylabel('magnitude [dB]'), title(['Spectrum (DFT) of time signal ''' self.label '''']);
					axis tight;
                    subplot(212);
                    stem(freq, 180/pi*unwrap(angle(val)));
                    set(gca,'xscal','log');
                    ylabel('phase [�]'), xlabel('f  [Hz]');
                    axis tight;
                otherwise
                    error('I ended up in an undefined state. Shouldn''t happen.');
            end
        end
        
        function sig = mtimes(self,other)
            if isnumeric(self)
                m = self; sig = other;
            else
                m = other; sig = self;
            end
            sig.signal_.Value = sig.signal_.Value*m;
        end
    end
    
end

