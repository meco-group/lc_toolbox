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

classdef Weight
% Static class which groups all functions to construct convenient weighting
% functions for H-infinity/H2 synthesis.
%
% For the time being only LTI weighting functions are implemented,
% parameter-varying weights are planned. 
    
    methods (Static)
        function W = DC(dc,varargin)
            % Constructs a DC weighting function (a constant). 
            %
            % Parameters:
            %  dc : the DC value of inverse of the weight @type double
            %  varargin : may contain 
            %   - \c ''Ts'' or \c ''fs'' followed by a double specifying the sampling time 
            %   - char specifying how \c dc should be interpreted;
            %  either \c ''linear'' or \c ''db'' (default)
            % 
            % Return values:
            %  W : standard MATLAB transfer function defining the weighting
            %  function @type tf
            
            [dc_mode,~,Ts,~] = Weight.parseVarargin(varargin{:});
            
            switch(dc_mode)
                case 'linear'
                    % do nothing
                case 'db'
                    dc = 10^(dc/20);
            end

            dc = 1/dc; % Change the dc gain to the inverse  
            W = zpk([],[],dc);
            if ~(isempty(Ts) || Ts==0)
                W = c2d(W,Ts,'Tustin');
            end
        end
        
        function W = RO(wco,order,rolloff,varargin)
            % Constructs a weight to enforce roll-off. Not 
            % interesting for the end user to use - use \c Weight.LF or 
            % \c Weight.HF instead.  
            %
            % Parameters:
            %  wco : the desired cross-over frequency @type double
            %  order : the order/steepness or the roll-off @type double
            %  rolloff : char specifying whether you'd like high-frequency
            %  or low-frequency roll-off; use \c ''HF'' or \c ''LF''
            %  respectively @type char
            %  varargin : may contain 
            %   - \c ''Ts'' or \c ''fs'' followed by a double specifying the sampling time 
            %   - char specifying how \c wco should be interpreted;
            %  either \c ''rad/s'' or \c ''Hz'' (default)
            %   - double specifying the gain at low-frequencies (if
            %   \c rolloff is \c ''LF'') or at high frequencies (if \c
            %   rolloff is \c ''HF'')
            %   - char specifying how the latter value should be interpreted;
            %  either \c ''linear'' or \c ''db'' (default)
            % 
            % Return values:
            %  W : standard MATLAB transfer function defining the weighting
            %  function @type tf
            
            [dc_mode,wco_mode,Ts,remainder] = Weight.parseVarargin(varargin{:});
            par = wco;
            switch(wco_mode)
                case 'Hz'
                    wco = wco*2*pi;
                case 'rad/s'
                    % do nothing
            end
            
            numeric = cellfun(@(x) isnumeric(x),remainder);
            switch(sum(numeric))
                case 0
                    % Make unstable or improper roll off weight
                    B = wco^order;
                    A = zeros(1,order+1); A(1) = 1;
                    switch(rolloff)
                        case 'LF'
                            switch isnumeric(wco)
                                case 0
                                    W = TFmod(B,A,par,Ts);
                                case 1
                                    W = tf(B,A);
                            end
                        case 'HF'
                            switch isnumeric(wco)
                                case 0
                                    W = TFmod(A,B,par,Ts);
                                case 1
                                    W = tf(A,B);
                            end
                    end

                case 1
                    % Make weight with flatgain
                    flatgain = remainder{numeric};
                    switch(dc_mode)
                        case 'linear'
                            % do nothing
                        case 'db'
                            flatgain = 10^(flatgain/20);
                    end
                    flatgain = 1/flatgain;

                    % Construct weight
                    switch(rolloff)
                        case 'LF'
                            w0 = wco/((flatgain^2-1)^(1/(2*order)));
                            fpass = 'low';
                        case 'HF'
                            w0 = wco*((flatgain^2-1)^(1/(2*order)));
                            fpass = 'high';
                    end

                    if w0<0
                        error('The break frequency of the FO filter is negative. This is probabely because your (linear) dc gain is below 1.')
                    end  
                    
                    [b,a] = butter2(order,w0,fpass,'s');
                    switch isnumeric(wco)
                        case 0 
                            W = flatgain*TFmod(b,a,par,Ts);
                        case 1
                            W = flatgain*tf(b,a);
                    end
                otherwise
                    error('RO cannot handle more than 1 numeric input');
            end
            
            if ~(isempty(Ts) || Ts==0)
                W = c2d(W,Ts,'Tustin');
            end
        end
        
        function [W] = LF(fco,order,varargin)
            % Constructs a weight to enforce low-frequency roll-off.
            %
            % Parameters:
            %  fco : the desired cross-over frequency @type double
            %  order : the order/steepness or the roll-off @type double
            %  varargin : may contain 
            %   - \c ''Ts'' or \c ''fs'' followed by a double specifying the sampling time 
            %   - char specifying how \c fco should be interpreted;
            %  either \c ''rad/s'' or \c ''Hz'' (default)
            %   - double specifying the gain at low frequencies, \c 0 if
            %   omitted (unstable weight)
            %   - char specifying how the latter value should be interpreted;
            %  either \c ''linear'' or \c ''db'' (default)
            % 
            % Return values:
            %  W : standard MATLAB transfer function defining the weighting
            %  function @type tf
            if(nargin>2)
                W = Weight.RO(fco,order,'LF',varargin{:});
            else
                W = Weight.RO(fco,order,'LF',varargin{:});
            end
        end
        
        function [W] = HF(fco,order,varargin)
            % Constructs a weight to enforce high-frequency roll-off.
            %
            % Parameters:
            %  fco : the desired cross-over frequency @type double
            %  order : the order/steepness or the roll-off @type double
            %  varargin : may contain 
            %   - \c ''Ts'' or \c ''fs'' followed by a double specifying the sampling time 
            %   - char specifying how \c fco should be interpreted;
            %  either \c ''rad/s'' or \c ''Hz'' (default)
            %   - double specifying the gain at high frequencies, \c 0 if
            %   omitted (improper weight)
            %   - char specifying how the latter value should be interpreted;
            %  either \c ''linear'' or \c ''db'' (default)
            % 
            % Return values:
            %  W : standard MATLAB transfer function defining the weighting
            %  function @type tf
            if(nargin>2)
                W = Weight.RO(fco,order,'HF',varargin{1},'db','Hz');
            else
                W = Weight.RO(fco,order,'HF','db','Hz');
            end
        end
        
        function [W] = PV(fco,order,varargin)
            % Constructs a weight to enforce parameter varying weights.
            %
            % Parameters:
            %  fco : the desired cross-over frequency @type double
            %  order : the order/steepness or the roll-off @type double
            %  varargin : may contain 
            %   - \c ''Ts'' or \c ''fs'' followed by a double specifying the sampling time 
            %   - char specifying how \c fco should be interpreted;
            %  either \c ''rad/s'' or \c ''Hz'' (default)
            %   - double specifying the gain at high frequencies, \c 0 if
            %   omitted (improper weight)
            %   - char specifying how the latter value should be interpreted;
            %  either \c ''linear'' or \c ''db'' (default)
            % 
            % Return values:
            %  W : standard MATLAB transfer function defining the weighting
            %  function @type tf
            if(nargin>2)
                W = Weight.RO(fco,order,'HF',varargin{1},'db','Hz');
            else
                W = Weight.RO(fco,order,'HF','db','Hz');
            end
        end
        
        function [dc_mode,wco_mode,Ts,remainder] = parseVarargin(varargin)
        % Parses input arguments and sets default values for the 
        % construction of weighting functions. Not interesting for the end 
        % user to use. 
        %
        % Parameters:
        %  varargin : may contain 
        %   - \c ''Ts'' or \c ''fs'' followed by a double specifying the sampling time 
        %   - char specifying how frequency values should be interpreted;
        %  either \c ''rad/s'' or \c ''Hz'' (default)
        %   - char specifying how gains should be interpreted;
        %  either \c ''linear'' or \c ''db'' (default)
        %   - other values are put in cell \c remainder 
        % 
        % Return values:
        %  dc_mode : char specifying how gains should be interpreted @type
        %  char
        %  wco_mode : char specifying how frequencies should be interpreted
        %  @type char
        %  Ts : sampling time @type double
        %  remainder : cell containing information that could not be casted
        %  into one of the former properties @type cell
        %  parametric : contains the fact if the weight is going to be
        %  parameteric or not
        
            %Set defaults
            dc_mode = 'db';
            wco_mode = 'Hz';
            Ts = 0;
            
            %Do stringwise checks
            checkLinear = cellfun(@(x) strcmp(x,'linear'),varargin);
            checkRads = cellfun(@(x) strcmp(x,'rad/s'),varargin);
            checkTs = cellfun(@(x) strcmp(x,'Ts'),varargin);
            checkfs = cellfun(@(x) strcmp(x,'fs'),varargin);
            
            % Parse checks
            if(any(checkLinear))
                dc_mode = 'linear';
            end
            if(any(checkRads))
                wco_mode = 'rad/s';
            end
            
            % Check if the weight is supposed to be discrete time
            if any(checkTs)
                if(find(checkTs) == length(checkTs))
                    error('A discrete weight requires ''Ts'' and a specified sample time');
                else
                    Ts = varargin{circshift(checkTs,[0 1])};
                    if(~isnumeric(Ts))
                        error('A discrete weight requires ''Ts'' and a specified sample time');
                    end
                end
            elseif any(checkfs)
                if(find(checkfs) == length(checkfs))
                    error('A discrete weight requires ''fs'' and a specified sample frequency in [Hz] or [rad/s] depending on extra option');
                else
                    fs = varargin{circshift(checkfs,[0 1])};
                    if(~isnumeric(fs))
                        error('A discrete weight requires ''fs'' and a specified sample frequency in [Hz] or [rad/s] depending on extra option');
                    end
                    Ts = 1/fs;
                    if(strcmp(wco_mode,'rad/s'))
                        Ts = Ts/(2*pi);
                    end
                end
            end
            
            % Remove part of the arguments to get the remainder
            varargin(checkTs | circshift(checkTs,[0 1]) | ...
                     checkfs | circshift(checkfs,[0 1]) |...
                     checkLinear | checkRads) = [];
            remainder = varargin;
        end
    end
end

