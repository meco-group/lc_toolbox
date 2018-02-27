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

classdef(InferiorClasses = {}) Norm
    %NORM
    %   Supports the definition of a multivariate signal norm. A norm
    %   consists of a weighting function W, an abstract transferfunction TF
    %   indicating the input-output relation to be weighted. W_Place
    %   indicates if the weight is an input(1) or output(0) weight.
    %   When dealing with a multivariable norm, W, W_place and TF are
    %   cellarrays of equal size.
    %   type specifies the type of the norm and can either be 'inf' or '2'.
    %   multiplier specifies an extra multiplier for the norm and can be
    %   used to weight the norm (e.g. in a constraint) or readjust the
    %   magnitude after calculating an optimal 'gamma' after optimization
    %   of the objective.
    %   Concatenation of 2 norms in the vertical direction results in a new
    %   norm. Horizontal concatenation of 2 norms results in an array of
    %   separate norms. Adding norms is allowed only when in an objective.
    %   Norm <= alpha results in a weighting of the norm with a factor
    %   alpha and is only allowed when defining a constraint.
    
    properties
        p           = Inf   % Type of the norm (2 or Inf)
        W_in        = [];   % Input weight of the channel ch_p
        W_out       = [];   % Output weight for the channel ch_p
        ch_p        = [];   % Input and output of the channel
        wscale      = [];   % Scaling of the weight in the optimization objective 
    end
    
    methods
        function self = Norm(varargin)
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
                        self.W_out = eye(length(self.ch_p.out));
                    elseif isa(varargin{2},'Channel')
                        self.ch_p = varargin{2};
                        self.W_in = eye(length(self.ch_p.in));
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
        
        function b = isHinf(self)
            b = (self.p == Inf);
        end
        
        function b = isH2(self)
            b = (self.p == 2);
        end
        
        function b = isoutput(self)
            b = arrayfun(@(x) ~isnumeric(x.W_out),self);
        end
        
        function b = isinput(self)
            b = arrayfun(@(x) ~isnumeric(x.W_in),self);
        end
        
        function b = isstable(self)
            b = (~isoutput(self) || isstable(self.W_out)) && (~isinput(self) || isstable(self.W_in));
%             b = isstable(self.W_out) && isstable(self.W_in);
        end
                
        function b = isparametric(self)
            b = false;
            if ~isnumeric(self.W_in), b = b | isparametric(self.W_in); end
            if ~isnumeric(self.W_out), b = b | isparametric(self.W_out);end
        end
        
        function s = scale(self)
            if length(self) == 1
                s = 1;
                if ~isempty(self.wscale)
                    s = self.wscale;
                end
            else
                s = arrayfun(@scale,self);
            end
        end
        
        function n = le(n,v)
            n = dealscale(n);
            n = n*(1/v);
        end
        
        function n = mtimes(a,b)
            if(isa(a,'Norm'))
                n = a; v = b; pre = true;
            else
                n = b; v = a; pre = false;
            end
            if isnumeric(v)
                if isempty(n.wscale)
                    n.wscale = v;
                else
                    n.wscale = n.wscale*v;
                end
            else
                v = fromstd(v);
                if pre
                    n.W_in = v;
                else
                    n.W_out = v;
                end
            end
        end
        
        function n = times(a,b)
            n = arrayfun(@(x,y)mtimes(x,y),a,b,'un',0);
            n = horzcat(n{:});
        end
        
        function self = dealscale(self,scale)
            if length(self) == 1
                if nargin < 2
                    scale = self.wscale;
                    self.wscale = [];
                end
                assert(any(length(scale)==[0,1]), 'Scale must be of the same size as norm');
            
                if ~isempty(scale)
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
                end
            else
                if nargin < 2
                    for k = 1:length(self)
                        self(k) = dealscale(self(k));
                    end
                else
                    for k = 1:length(self)
                        self(k) = dealscale(self(k),scale(k));
                    end
                end
            end
        end
                
        %% Norm stacking function
        function self = vertcat(varargin)
            self = varargin{1};
            for k = 2:nargin
                self = stack(self,varargin{k});
            end
        end
        
        function self = stack(self,other)
            % Stack to norms together
            if self.p ~= other.p
                error('cannot concatenate 2 norms of different types.');
            else % norms can be concatenated
                ch_p_ = Channel([self.ch_p.in;other.ch_p.in],[self.ch_p.out;other.ch_p.out]);
                W_out_ = blkdiag(self.W_out,other.W_out);
                W_in_ = blkdiag(self.W_in,other.W_in);
                self = Norm(ch_p_);
                self.W_out = W_out_;
                self.W_in = W_in_;
                self = setP(self,other.p);
            end
        end
        
        function self = norm(self,p)
            if nargin > 1
                self = setP(self,p);
            end
        end
        
        function [norms,performances] = unstack(self,performance)
            if all(self.isoutput())
                [norms,performances] = unstackout(self,performance);
            else
                [norms,performances] = unstackin(self,performance);
            end
        end
        
        function narray = plus(n1,n2)
            narray = [n1,n2];
        end
        
        function [w,b] = weight(self,out,in)
            assert(length(out)==1 && length(in)==1,'Only implemented for siso channels');
            if length(self) == 1
                self = combine(self);
                w = {};
                b = false;
                [lio,loco] = ismember(self.ch_p.out,out);
                [lii,loci] = ismember(self.ch_p.in,in);
                if all(lii) && all(lio)
                    if isinput(self) && isoutput(self)
                        assert(all([size(wout),size(win)]==1),'This function only works for siso');
                        w = {self.W_out(:,loco)*self.W_in(loci,:)};
                        b = true;
                    elseif isinput(self)
                        w = {self.W_in(loci,:)};
                        b = true;
                    elseif isoutput(self)
                        w = {self.W_out(:,loco)};
                        b = true;
                    end
                end
            else
                [w,b] = arrayfun(@(x)weight(x,out,in),self,'un',0);
                w = horzcat(w{:}); b = any(horzcat(b{:}));
            end
        end
        
        function [w,b] = invweight(self,out,in)
            [w,b] = weight(self,out,in);
            if b
                w = cellfun(@inv,w,'un',0);
            end
        end
        
        function [w,b] = sigma(self,out,in,varargin)
            [w,b] = weight(self,out,in); f = {};
            if b
                [sv,f] = cellfun(@(x)sigma(x,varargin{:}),w,'un',0);
                w = cellfun(@(x,y)FRDmod(x(1,:),y),sv,f,'un',0);
            end
        end
        
        function [w,b] = invsigma(self,out,in,varargin)
            [w,b] = sigma(self,out,in,varargin{:});
            if b
                w = cellfun(@inv,w,'un',0);
            end
        end
        
        function self = combine(self)
            function [c,S] = getS(a)
                [c,~,ic] = unique(a);
                S = double(cell2mat(arrayfun(@(x)eq(x,ic),1:length(c),'un',0)))';
            end
            if numel(self) == 1
                [self.ch_p.in,Si] = getS(self.ch_p.in);
                if isinput(self)
                    self.W_in = Si*self.W_in;
                else
                    self.W_in = eye(length(self.ch_p.in));
                end
                [self.ch_p.out,So] = getS(self.ch_p.out);
                if isoutput(self) 
                    self.W_out = self.W_out*transpose(So);
                else
                    self.W_out = eye(length(self.ch_p.out));
                end
            else
                c = arrayfun(@combine,self,'un',0);
                self = horzcat(c{:});
            end
        end
    end
        
    methods(Access=private)        
        function self = setP(self,p)
            if (p==2)||(p==inf)
                self.p = p;
            else
                error('System norm should be 2 or inf');
            end
        end
        
                
%         function self = stackout(self,other)
%             newout = setdiff(other.TF.out,self.TF.out);
%             self.TF.out = [self.TF.out;newout];
% 
%             self.W = blkdiag(self.W,ss(zeros(length(other.TF.in),length(newout))));
%             [~,i] = ismember(other.TF.out,self.TF.out);
%             self.W((length(self.TF.in)+1):end,i) = other.W;
%             
%             self.TF.in = [self.TF.in;other.TF.in];
%         end
%         
%         function self = stackin(self,other)
%             self = transpose(stackout(transpose(self),transpose(other)));
%         end
        
        function [norms,performances] = unstackout(self, performance)
            if size(self.W,2)>1
                self = self.unify();
                active = entrynorm(self.W) > eps;
                N = size(active,2);

                % allocate memory
                norms(1,N) = Norm();
                performances = cell(N,1);

                for k = 1:N % iterate over columns
                    % Save part of the weight as new norm
                    norms(1,k) = Norm(Channel(self.TF.in(active(:,k)),self.TF.out(k)), self.W(active(:,k),k), self.position, self.p);
                    if ~isempty(performance)
                        performances{k,1} = performance(k,:);
                    end
                end
            else
                norms = self;
                performances = {performance};
            end
        end
        
        function [norms,performances] = unstackin(self, performance)
            [norms,performances] = unstackout(transpose(self),performance);
            norms = transpose(norms);
        end
        
        function self = transpose(self)
            if length(self)>1
                for k=1:length(self)
                    self(k) = transpose(self);
                end
            else
                t = self.TF.in;
                self.TF.in = self.TF.out;
                self.TF.out = t;

                self.W = transpose(self.W);
            end
        end
    end
end

