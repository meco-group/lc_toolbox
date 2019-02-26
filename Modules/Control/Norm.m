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

classdef Norm < Specification
    % A Norm defines the weighted norm (@f$\mathcal{H}_2@f$ or @f$\mathcal{H}_\infty@f$) of a certain performance
    % channel as an objective, or as a constraint (becomes NormConstraint) 
    % when an upper bound is imposed.
    
    properties
        p           = Inf   % type of the norm: 2 for @f$\mathcal{H}_2@f$ norm, Inf for @f$\mathcal{H}_\infty@f$ norm (default)
        W_in        = [];   % input weight of the channel ch_p
        W_out       = [];   % output weight for the channel ch_p
        ch_p        = [];   % input and output of the channel
        wscale      = [];   % scaling of the weight in the optimization objective 
    end
    
    methods
        
        function self = Norm(varargin)
        % Constructor for Norm objects.
        % 
        % Parameters:
        %   varargin : should have on of these four structures:
        %   -# an output weight (\c numlti or Model), followed by the performance
        %      channel (\c Channel) , and the input weight (\c numlti or Model) as last argument 
        %   -# an output weight (\c numtli or Model), followed by the
        %      performance channel (\c Channel)
        %   -# the performance channel (\c Channel), followed by an input weight (\c numlti or Model)
        %   -# a performance channel (\c Channel) (no weighting)
        %
        % Return values:
        %  self : the weighted norm @type Norm
            
            if nargin > 3
                self.p = varargin{end};
                varargin(end) = [];
            end
            switch length(varargin)
                case 3 % Wout - channel - Win
                    self.ch_p = varargin{2};
                    self.W_in = fromstd(varargin{3});
                    self.W_out = fromstd(varargin{1});
                case 2 
                    if isa(varargin{1},'Channel') 
                        self.ch_p = varargin{1};
                        self.W_in = fromstd(varargin{2});
                        self.W_out = fromstd(eye(length(self.ch_p.out)));
                    elseif isa(varargin{2},'Channel')
                        self.ch_p = varargin{2};
                        self.W_in = fromstd(eye(length(self.ch_p.in)));
                        self.W_out = fromstd(varargin{1});
                    elseif isa(varargin{1},'Norm') && isnumeric(varargin{2})
                        self = varargin{1};
                        self.p = varargin{2};
                    else
                        error('At least one of the input arguments should be a channel');
                    end
                case 1
                    assert(isa(varargin{1},'Channel'),'The argument should be a channel');
                    self.ch_p = varargin{1};
                    self.W_out = eye(length(self.ch_p.out));
                    self.W_in = eye(length(self.ch_p.in));
            end
            
            assert((size(self.W_in,1)==length(self.ch_p.in)),'Input weight dimensions do not correspond to the channel');
            assert((size(self.W_out,2)==length(self.ch_p.out)),'Output weight dimensions do not correspond to the channel');
            
        end
        
        function disp(self)
        % Print information
            if isH2(self)
                type = 'H2';
            else
                type = 'H-infinity';
            end
            inputs = ['[',strjoin(arrayfun(@(x)num2str(x.UUID),self.ch_p.in','un',0),';'),']'];
            outputs = ['[',strjoin(arrayfun(@(x)num2str(x.UUID),self.ch_p.out','un',0),';'),']'];
            
            disp([type, ' norm on channel ', self.ch_p.name, ' : ' inputs, ' ->  ', outputs]);
        end
        
        function b = isHinf(self)
        % Checks whether the norm is an @f$\mathcal{H}_\infty@f$ norm.
        %
        % Return values:
        %  b : true if the norm is an @f$\mathcal{H}_\infty@f$ norm @type logical
            b = (self.p == Inf);
        end
        
        function b = isH2(self)
        % Checks whether the norm is an @f$\mathcal{H}_2@f$ norm.
        %
        % Return values:
        %  b : true if the norm is an @f$\mathcal{H}_2@f$ norm @type logical
            b = (self.p == 2);
        end
        
        function b = isoutput(self)
        % Checks whether the norm has a non-static output weight.
        %
        % Return values:
        %  b : true if the norm has a non-static output weight @type logical
            b = arrayfun(@(x) ~isequal(x.W_out,SSmod(eye(size(x.W_out,2)))),self);
        end
        
        function b = isinput(self)
        % Checks whether the norm has a non-static input weight.
        %
        % Return values:
        %  b : true if the norm has a non-static input weight @type logical
            b = arrayfun(@(x) ~isequal(x.W_in,SSmod(eye(size(x.W_in,1)))),self);
        end
        
        function b = isstable(self)
        % Checks whether the weights are stable.
        %
        % Return values:
        %  b : true if both input and output weights are stable @type logical
            b = (~isoutput(self) || isstable(self.W_out)) && (~isinput(self) || isstable(self.W_in));
        end
                
        function b = isparametric(self)
        % Checks whether the weights are parameter-dependent. 
        %
        % Return values:
        %  b : true if the weights are parameter-dependent @type logical
            b = false;
            if ~isnumeric(self.W_in), b = b | isparametric(self.W_in); end
            if ~isnumeric(self.W_out), b = b | isparametric(self.W_out);end
        end
        
        function b = compin(self,other)
        % Checks whether the two norms have the same input and same input weight
        %
        % Return values:
        %  b : true if the norms share the same input and weight @type logical
            bw = sscomp(self.W_in, other.W_in);
            bs = isalias(self.ch_p.in,other.ch_p.in);
            b = bw && bs;
        end
        
        function b = compout(self,other)
        % Checks whether the two norms have the same output and same output weight
        %
        % Return values:
        %  b : true if the norms share the same output and weight @type logical
            bw = sscomp(self.W_out, other.W_out);
            bs = isalias(self.ch_p.out,other.ch_p.out);
            b = bw && bs;
        end
        
        function s = scale(self)
        % Returns the scale of the norm, e.g. @f$\alpha@f$ for
        % @f$\alpha \|\cdot\|_{2/\infty}@f$.
        %
        % Return values:
        %  s : the scale factor of the norm @type double
            if length(self) == 1
                s = 1;
                if ~isempty(self.wscale)
                    s = self.wscale;
                end
            else
                s = arrayfun(@scale,self);
            end
        end
        
        function n = lt(n,v)
        % Creates a NormConstraint based on a Norm by the operator <. 
        % Overloads MATLAB's \c < operator for Norm objects. 
        %
        % Parameters:
        %   n : the norm that is constrained @type Norm
        %   v : the upper bound on the norm @type double
        %
        % Return values:
        %  s : the constraint @type NormConstraint
            n = le(n,v);
        end
        
        function nc = le(n,v)
        % Creates a NormConstraint based on a Norm by the operator <=. 
        % Overloads MATLAB's \c <= operator for Norm objects. 
        %
        % Parameters:
        %   n : the norm that is constrained @type Norm
        %   v : the upper bound on the norm @type double
        %
        % Return values:
        %  s : the constraint @type NormConstraint
            nc = NormConstraint(n,v);
        end
        
        function n = mtimes(a,b)
        % Multiplies a norm with a scalar or a weight. 
        % Overloads MATLAB's \c * for Norm objects. 
        %
        % Parameters:
        %   a : either a Norm, a weight (\c numlti or Model), or a scalar 
        %   b : either a Norm, a weight (\c numlti or Model), or a scalar 
        %
        % Return values:
        %  n : the scaled or weighted norm @type Norm
            if(isa(a,'Norm'))
                n = a; v = b; pre = true;
            else
                n = b; v = a; pre = false;
            end
            if isnumeric(v) && isscalar(v)
                if isempty(n.wscale)
                    n.wscale = v;
                else
                    n.wscale = n.wscale*v;
                end
            else
                v = fromstd(v);
                if pre
                    if isempty(n.W_in)
                        n.W_in = v;
                    else
                        n.W_in = n.W_in*v;
                    end
                else
                    if isempty(n.W_out)
                        n.W_out = v;
                    else
                        n.W_out = v*n.W_out;
                    end
                end
            end
        end
        
        function n = times(a,b)
        % Multiplies a norm with a scalar or a weight. 
        % Overloads MATLAB's \c .* for Norm objects. 
        %
        % Parameters:
        %   a : either a Norm, a weight (\c numlti or Model), or a scalar 
        %   b : either a Norm, a weight (\c numlti or Model), or a scalar 
        %
        % Return values:
        %  n : the scaled or weighted norm @type Norm
            n = arrayfun(@(x,y)mtimes(x,y),a,b,'un',0);
            n = horzcat(n{:});
        end
        
        function self = vertcat(varargin)
        % Stacks multiple norms with the same inputs into one norm. 
        % Overloads MATLAB's \c vertcat for Norm objects. 
        %
        % Parameters:
        %   varargin: set of Norm objects 
        %
        % Return values:
        %  n : the stacked norm @type Norm
            if all(cellfun(@(x) isa(x,'Norm'),varargin))
                self = varargin{1};
                for k = 2:nargin
                    self = stackvert(self,varargin{k});
                end
            else
                self = vertcat@Specification(varargin{:});
            end
        end
        
        function self = horzcat(varargin)
        % Stacks multiple norms with the same outputs into one norm. 
        % Overloads MATLAB's \c horzcat for Norm objects. 
        %
        % Parameters:
        %   varargin: set of Norm objects 
        %
        % Return values:
        %  n : the stacked norm @type Norm
            if all(cellfun(@(x) isa(x,'Norm'),varargin))
                self = varargin{1};
                for k = 2:nargin
                    self = stackhorz(self,varargin{k});
                end
            else
                self = horzcat@Specification(varargin{:});
            end
        end
        
        function self = blkdiag(varargin)
        % Stacks multiple norms into one norm. 
        % Overloads MATLAB's \c blkdiag for Norm objects. 
        %
        % Parameters:
        %   varargin: set of Norm objects 
        %
        % Return values:
        %  n : the stacked norm @type Norm
            if all(cellfun(@(x) isa(x,'Norm'),varargin))
                self = varargin{1};
                for k = 2:nargin
                    self = stack(self,varargin{k});
                end
            else
                error('blkdiag for norms is only applicable if all input arguments are of type Norm.');
            end
        end
        
        function self = stackvert(self,other)
        % Stacks multiple norms with the same inputs into one norm. Related to Channel::vertcat.
        %
        % Parameters:
        %   self: the first norm @type Norm
        %   other: the second norm @type Norm
        %
        % Return values:
        %  self : the stacked norm @type Norm
            if isempty(self)
                self = other;
            elseif isempty(other)
                return; 
            else
                assert(sscomp(self.W_in,other.W_in), 'If you want to concatenate norms vertically, their input weights should be equal. Use blkdiag instead.');
                if self.p ~= other.p
                    error('Cannot stack 2 norms of different types.');
                else
                    ch_p_ = [self.ch_p;other.ch_p];
                    W_out_ = blkdiag(self.W_out,other.W_out);
                    W_in_ = self.W_in;
                    self = Norm(ch_p_);
                    self.W_out = W_out_;
                    self.W_in = W_in_;
                    self = setP(self,other.p);
                end
            end
        end
        
        function self = stackhorz(self,other)
        % Stacks multiple norms with the same outputs into one norm. Related to Channel::horzcat.
        %
        % Parameters:
        %   self: the first norm @type Norm
        %   other: the second norm @type Norm
        %
        % Return values:
        %  self : the stacked norm @type Norm
            if isempty(self)
                self = other;
            elseif isempty(other)
                return; 
            else
                assert(sscomp(self.W_out,other.W_out), 'If you want to concatenate norms horizontally, their output weights should be equal. Use blkdiag instead.');
                if self.p ~= other.p
                    error('Cannot stack 2 norms of different types.');
                else 
                    ch_p_ = [self.ch_p other.ch_p];
                    W_out_ = self.W_out;
                    W_in_ = blkdiag(self.W_in,other.W_in);
                    self = Norm(ch_p_);
                    self.W_out = W_out_;
                    self.W_in = W_in_;
                    self = setP(self,other.p);
                end
            end
        end
        
        function self = stack(self,other)
        % Stacks multiple norms into one norm. Related to Channel::blkdiag.
        %
        % Parameters:
        %   self: the first norm @type Norm
        %   other: the second norm @type Norm
        %
        % Return values:
        %  self : the stacked norm @type Norm
            if isempty(self)
                self = other;
            elseif isempty(other)
                return; 
            else
                if self.p ~= other.p
                    error('Cannot stack 2 norms of different types.');
                else 
                    ch_p_ = blkdiag(self.ch_p, other.ch_p);
                    W_out_ = blkdiag(self.W_out,other.W_out);
                    W_in_ = blkdiag(self.W_in,other.W_in);
                    self = Norm(ch_p_);
                    self.W_out = W_out_;
                    self.W_in = W_in_;
                    self = setP(self,other.p);
                end
            end
        end
        
        function self = setnorm(self,p)
        % Sets the norm type. 
        %
        % Parameters:
        %   self: the norm @type Norm
        %   p: the norm type (2 or Inf) @type double
        %
        % Return values:
        %  self : the resulting norm @type Norm
            if nargin > 1
                self = setP(self,p);
            end
        end
        
        function self = getnorm(self)
        % Returns itself (used for the sake of simplicity). 
        %
        % Parameters:
        %   self: the norm @type Norm
        %
        % Return values:
        %  self : the norm @type Norm
        end
        
        function self = setWout(self,W_out)
        % Sets the output weight of the norm.
        %
        % Parameters:
        %   self: the norm @type Norm
        %   W_out: the output weight of the norm @type Model
        %
        % Return values:
        %  self : the resulting norm @type Norm
            self.W_out = W_out; 
        end
        
        function self = setWin(self,W_in)
        % Sets the input weight of the norm.
        %
        % Parameters:
        %   self: the norm @type Norm
        %   W_in: the input weight of the norm @type Model
        %
        % Return values:
        %  self : the resulting norm @type Norm
            self.W_in = W_in; 
        end
        
        function w = getWin(self)
        % Returns the input weight
        %
        % Parameters:
        %   self: the norm @type Norm
        %
        % Return values:
        %  w : the input weight @type DynamicSystem
        
            w = self.W_in;
        end
        
        function w = getWout(self)
        % Returns the output weight
        %
        % Parameters:
        %   self: the norm @type Norm
        %
        % Return values:
        %  w : the output weight @type DynamicSystem
        
            w = self.W_out;
        end
        
        function self = setin(self,signal)
        % Sets the input signals of the norm channel.
        %
        % Parameters:
        %   self: the norm @type Norm
        %   signal: the input signals of the norm channel @type Signal
        %
        % Return values:
        %  self : the resulting norm @type Norm
            self.ch_p.in = signal;
        end
        
        function self = setout(self,signal)
        % Sets the output signals of the norm channel.
        %
        % Parameters:
        %   self: the norm @type Norm
        %   signal: the output signals of the norm channel @type Signal
        %
        % Return values:
        %  self : the resulting norm @type Norm
            self.ch_p.out = signal;
        end
        
        function self = setchannel(self,ch)
        % Sets the channel of the norm.
        %
        % Parameters:
        %   self: the norm @type Norm
        %   ch: the channel of the norm @type Channel
        %
        % Return values:
        %  self : the resulting norm @type Norm
            self.ch_p = ch;
        end
        
        function ch = getchannel(self)
        % Returns the channel of the norm.
        %
        % Parameters:
        %   self: the norm @type Norm
        %
        % Return values:
        %  ch : the channel of the norm @type Channel
            ch = self.ch_p;
        end
        
        function p = normtype(self)
        % Returns the norm type (2 or Inf).
        %
        % Parameters:
        %   self: the norm @type Norm
        %
        % Return values:
        %  self : the norm type (2 or Inf)
            p = self.p;
        end
        
        function sum = plus(varargin)
        % Returns a set of Norm objects as a cell.
        %
        % Parameters:
        %  varargin : can contain all Norm objects and cells
        %  containing a set of Norm objects 
        %
        % Return values:
        %  c : cell containing all norms @type cell
            isnorm = cellfun(@(x) isa(x,'Norm'), varargin);
            iscellofnorms = cellfun(@(x) isa(x,'cell') && all(cellfun(@(y) isa(y,'Norm'),x)), varargin);
            
            if any(~(isnorm | iscellofnorms)) 
                sum = plus@Specification(varargin{:}); 
            else
                norms = varargin(isnorm);
                if ~isempty(varargin(iscellofnorms))
                    unpacked_norms = vertcat(varargin{iscellofnorms});
                else
                    unpacked_norms = cell(1,0);
                end
            
                sum = [norms(:) ; unpacked_norms(:)];
            end
        end
    end
        
    methods(Access=private)        
        function self = setP(self,p)
        % Sets the norm type. 
        %
        % Parameters:
        %   self: the norm @type Norm
        %   p: the norm type (2 or Inf) @type double
        %
        % Return values:
        %  self : the resulting norm @type Norm
            if (p==2)||(p==inf)
                self.p = p;
            else
                error('System norm should be 2 or Inf.');
            end
        end
    end
    
    methods(Static)
        function self = dealscale(self,scale,io)
        % Includes the scale of the norm in the weights and resets the
        % norm's scale to 1 (or empty). Is a static method since it is also
        % possible to call it on cells. 
        %
        % Parameters:
        %   self: the norm (\c Norm), or a cell of norms (\c cell)
        %   scale: optional scale, which the current scale is multiplied
        %   with; 1 if omitted (\c double or \c cell)
        %
        % Return values:
        %  self : the same (set of) norm(s), but with the scale included in
        %  the weights (\c Norm or \c cell)

            if nargin < 2
                scale = num2cell(ones(size(self)));
            else
                assert(isa(scale,'cell') && length(scale)==length(self), 'scale should be a cell with as many scales as norms.');
            end

            if isa(self,'cell')
                self = cellfun(@(x,y) {Norm.dealscale(x,{y})}, self, scale(:));
            else
                assert(isa(self,'Norm'), 'dealscale can only operate on Norm objects.');
                if ~isempty(self.wscale)
                    scale = self.wscale*scale{:};
                else
                    scale = scale{:};
                end
                if isinput(self) && isoutput(self)
                    self.W_in = self.W_in*sqrt(scale);
                    self.W_out = self.W_out*sqrt(scale);
                elseif isinput(self)
                    self.W_in = self.W_in*scale;
                elseif isoutput(self)
                    self.W_out = self.W_out*scale;
                else
                    self.W_out = fromstd(eye(length(self.ch_p.out))*scale);
                end 
                self.wscale = [];
            end
        end
        
        function [unorms,ui,i] = unique(norms,io)
        % Compute the unique weighted inputs or outputs of a series of norms 
        %
        % Parameters:
        %   norms: cell of norms (\c Norm or \c cell)
        %   io: string to set comparison: input ('in') or output ('out')
        %   (\c strings)
        %
        % Return values:
        %  unorms: cell of the unique norms @type cell
        %  i: index for each element in norms to the corresponding element 
        %  in unorms @type double
        
            % input data formatting
            switch(io)
                case 'in'
                    signals = cellfun(@(x) x.ch_p.in,norms,'un',0);
                    weights = cellfun(@(x) x.W_in,norms,'un',0);
                case 'out'
                    signals = cellfun(@(x) x.ch_p.out,norms,'un',0);
                    weights = cellfun(@(x) x.W_out,norms,'un',0);
            end

            % loop through norms to look for unique channels
            ui = 1; % indices of unique norms: norms(ui) == unorms
            i = ones(1,length(norms)); % indices so that unorms(i(k)) == norms(k)
            for k = 2:length(norms)
                bs = cellfun(@(x) isalias(signals{k},x),signals(ui));
                bw = cellfun(@(x) sscomp(weights{k},x),weights(ui));
                j = find(bs & bw);
                assert(length(j) < 2, 'There was more than one match for the norm which should not happen!');
                if ~isempty(j)
                    i(k) = j;
                else
                    ui(end+1) = k;
                    i(k) = length(ui);
                end
            end

            % construct output
            unorms = norms(ui);            
        end
        
        function [sys,connections,norms] = tofilter(norms,io)
            [~,ui,i] = Norm.unique(norms,io);
            
            switch(io)
                case 'in'
                    signals = cellfun(@(x) x.ch_p.in,norms(ui),'un',0);
                    weights = cellfun(@(x) transpose(x.W_in),norms(ui),'un',0);
                case 'out'
                    signals = cellfun(@(x) x.ch_p.out,norms(ui),'un',0);
                    weights = cellfun(@(x) x.W_out,norms(ui),'un',0);
            end
            
            filter = blkdiag(weights{:});
            signals = vertcat(signals{:});
            [signals,IA,IC] = unique(signals,'stable');
            C = zeros(length(IC),length(IA));
            ind = sub2ind(size(C),1:length(IC),IC');
            C(ind) = 1;
            filter = filter*C;
            if isnumeric(filter), filter = SSmod(filter); end
            
            switch(io)
                case 'in'
                    sys = IOSystem(transpose(filter));
                    connections = (sys.out == signals);
                    counter = 0;
                    for k = ui
                        norms{k} = setin(norms{k}, sys.in(counter+(1:length(norms{k}.ch_p.in))));
                        counter = counter + length(norms{k}.ch_p.in);
                    end
                    for k = 1:length(i), norms{k} = setin(norms{k},norms{ui(i(k))}.ch_p.in); end
                case 'out'
                    sys = IOSystem(filter);
                    connections = (sys.in == signals);
                    counter = 0;
                    for k = ui
                        norms{k} = setout(norms{k}, sys.out(counter+(1:length(norms{k}.ch_p.out))));
                        counter = counter + length(norms{k}.ch_p.out);
                    end
                    for k = 1:length(i), norms{k} = setout(norms{k},norms{ui(i(k))}.ch_p.out); end
            end
        end
    end
end

