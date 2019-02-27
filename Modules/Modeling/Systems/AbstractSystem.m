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

classdef AbstractSystem < handle
    %IOSystem Basic input-output system
    %   Most general representation of a system with input and output
    %   signals. The IOSystem contains other systems or/and models.
    
    properties
        in = []; % input signals
        out = []; % output signals
    end
    
    properties (Access=protected)
        content_ = {};% list of content of the IOSystem. Both IOSystem and models
    end
    
    properties (Access=public)        
        pin = [];
        pout = [];
        fbin = [];
        fbout = [];
        
        fbsigin = [];
        fbsigout = [];
        fb = [];
    end
    
    methods (Abstract)
        add(self,varargin)
        model(self,varargin)
    end
    
    methods (Access=protected)
        subsref_signal(self,in,out)
    end
    
    methods
        function self = AbstractSystem(varargin)
            if nargin > 0
                if isnumeric(varargin{1}) && isnumeric(varargin{2}) % Give number of inputs,outputs
                    self.in = Signal(varargin{1});
                    self.out = Signal(varargin{2});
                elseif isa(varargin{1},'Signal') && isa(varargin{2},'Signal')
                    self.in = varargin{1};
                    self.out = varargin{2};
                else
                    error('Unknown constructor');
                end
            end                
        end
        
        function c = content(self,index)
            if nargin > 1
                if length(index) > 1
                    c = self.content_(index);
                else
                    c = self.content_{index};
                end
            else
                c = self.content_;
            end
        end
        
        function b = isempty(self)
            b = isempty(self.content_);
        end
        
        function [self,restore] = empty(self)
            [self,restore] = subsystem_helper(self,[]);
        end
        
        function [self,restore] = remove(self,index)
            keep = (1:numod(self));
            keep(index) = [];
            [self,restore] = subsystem_helper(self,keep);
        end
        
        function self = replace(self,index,model)
            if ~iscell(model), model = {model}; end
            assert(all(cellfun(@(x) isa(x,'Model'),model)),'Can only add Model to a system');
            self.content_(index) = model;
        end
        
        function n = numod(self)
            n = length(self.content_);
        end
        
        function [varargout] = size(self,dim)
            siz = [length(self.out),length(self.in)]; %-length(self.Cin)

            if nargin==1 && nargout==1
                varargout{1} = siz;
            elseif nargin==1 && nargout==2
                varargout{1} = siz(1);
                varargout{2} = siz(2);
            elseif nargin==2 && nargout==1
                varargout{1} = siz(dim);
            elseif nargin==0
                disp(size(self,dim));
            else
                error('Unknown combination of input and output arguments.');
            end
        end
        
        function varargout = subsref(self,s)
            % SUBSREF reimplement subsref with signal indexing

            if strcmp(s(1).type,'()')
                if any(cellfun(@(x)isa(x,'Channel'),s(1).subs))
                    assert(length(s(1).subs)==1,'You can only give one channel as a subscript. If you want to retrieve multiple channels at once, stack them.');
                    ch = s(1).subs{:};
                    newsubs{1} = ch.out; 
                    newsubs{2} = ch.in;
                    s(1).subs = newsubs;
                    varargout = subsref(self,s);
                end
                if any(cellfun(@(x)isa(x,'Signal'),s(1).subs))
                    if isequal(s(1).subs{1},':');s(1).subs{1}=self.out;end
                    if isequal(s(1).subs{2},':');s(1).subs{2}=self.in;end
                    varargout = {subsref_signal(self,s(1).subs{2},s(1).subs{1})};
                end
                if length(s) > 1
                    [varargout{1:nargout}] = builtin('subsref',varargout{1},s(2:end));
                end
            else
                [varargout{1:nargout}] = builtin('subsref',self,s(:));
            end
        end
                        
        
        
        %% Utility
        function mod = std(self)
            self = model(self);
            mod = std(self.content_{1});
        end
        
        function disp(self)
            disptab(self,0);
        end
        
    end
    
    methods (Access=protected)
        function disptab(self,tabs)
            tabs = repmat('\t',[1,numel(tabs)]);
            
            fprintf([tabs,'System with %d inputs and %d outputs\n'],length(self.in),length(self.out));
            if ~isempty(self.in) && ~isempty(self.out)
                fprintf([tabs,'\tin: [%s], out: [%s]\n'],num2str([self.in.UUID]),num2str([self.out.UUID]));
            end
        end
    end
    
    methods(Access=public)
        function self = addfbsys(self,fbsys)
            assert(isa(fbsys,'IOSystem'),'Can only add IOSystem as hidden content');
            self.fbcontent = [self.fbcontent,{fbsys}];
        end
        
        function self = resolve(self,connections)
            % RESOLVE Resolve connections for iosystem
            % Make linear combination system
            for k = 1:numel(connections)
                if islincomb(connections{k})
                    extrain = transpose(horzcat(connections{k}.components));
                    self.fbsigin = [self.fbsigin;extrain];
                    self.fbsigout = [self.fbsigout;connections{k};extrain];
                    self.fb = blkdiag(self.fb,blkdiag(connections{k}.multipliers));
                    self.fb = [self.fb;[zeros(length(extrain),size(self.fb,2)-length(extrain)),eye(length(extrain))]];
                end
            end
            
            self.in = cellfun(@(x)x.in,self.content_,'UniformOutput',false);
            self.out = cellfun(@(x)x.out,self.content_,'UniformOutput',false);
            self.in = vertcat(self.in{:});
            self.out = vertcat(self.out{:});
            
            self = sepio(self,connections);
            
            [~,self.pin,self.fbout] = intersect(self.in,self.fbsigout);
            [~,self.pout,self.fbin] = intersect(self.out,self.fbsigin);
            self.pin(self.pin == 0) = [];
            self.pout(self.pout == 0) = [];
            self.fbin(self.fbin == 0) = [];
            self.fbout(self.fbout == 0) = [];
            
            self = addalias(self,connections);
        end
        
        function self = addalias(self,connections)
            % anonymous function
            function sigvec = addaliasvec(sigvec,conn)
                for c = 1:length(sigvec)
                    for d = 1:length(conn) 
                        if isalias(sigvec(c),conn{d})
                            sigvec(c) = addalias(sigvec(c),conn{d});
                            break;
                        end
                    end
                end
            end
            
            connections = Signal.group(connections);
            self.in = addaliasvec(self.in,connections);
            self.out = addaliasvec(self.out,connections);
            self.fbsigin = addaliasvec(self.fbsigin,connections);
            self.fbsigout = addaliasvec(self.fbsigout,connections);
        end
        
        function [self] = sepio(self,connections)
            connections = Signal.group(connections);
            
            for k = 1:length(connections)
                [membin,~] = ismember(self.in,connections{k});
                [membout,~] = ismember(self.out,connections{k});
                membin = find(membin);
                membout = find(membout);
                if (length(membin) == 1) && (length(membout) == 1)
                    self.fbsigin = [self.fbsigin;self.out(membout)];
                    self.fbsigout = [self.fbsigout;self.in(membin)];
                    self.fb = blkdiag(self.fb,1);
                elseif (length(membin) > 1) && (length(membout) <= 1)
                    self.fbsigin = [self.fbsigin;connections{k}];
                    self.fbsigout = [self.fbsigout;self.in(membin)];
                    self.fb = blkdiag(self.fb,ones(length(membin),1));
                elseif (length(membin) <= 1) && (length(membout) > 1)
                    self.fbsigin = [self.fbsigin;self.out(membout)];
                    self.fbsigout = [self.fbsigout;connections{k}];
                    self.fb = blkdiag(self.fb,ones(1,length(membout)));
                elseif isempty(membin) || isempty(membout)
                    % ignore this -> pure alias
                else
                    keyboard
                    error('InputOutput cluster discovered. Inform the developers.');
                end
            end
        end
        
        function [sys] = mtimes(s1,s2)
            assert(isa(s1,'AbstractSystem') && isa(s2,'AbstractSystem'), 's1 and s2 should be systems (not models)');
            assert(size(s2,1) == size(s1,2),'number of inputs does not match');
            conn = [s1.in == s2.out];
            sys = IOSystem(s1,s2,conn);
        end
    end
    
    methods (Access=public, Sealed=true)
        function varargout = bode(varargin)
            [varargin,labels] = stdargs(varargin,'freq');
                for i = 1:length(varargin)
                    if (numel(varargin{i}.SamplingGrid) == 1)
                    s =  fieldnames(varargin{i}.SamplingGrid);
                    param_grid{i} = varargin{i}.SamplingGrid.(s{1});
                        for j = 1:size(varargin{i},3)
                            [mag,phase,w] = bode(varargin{i}(:,:,j,:));
                            h1 = subplot(2,1,1);
                            hold(h1,'on')
                            h1.ColorOrderIndex = i;
                            plot3(w(:),param_grid{i}(j)*ones(size(w,1),1),mag(:));
                            h2 = subplot(2,1,2);
                            hold(h2,'on')
                            h2.ColorOrderIndex = i;
                            plot3(w(:),param_grid{i}(j)*ones(size(w,1),1),phase(:));
                        end
                    end
                end
                h1.View = [45,45];
                h2.View = [45,45];
                h1.XLabel.String = 'Frequency';
                h2.XLabel.String = 'Frequency';
                h1.ZLabel.String = 'Magnitude(dB)';
                h2.ZLabel.String = 'Phase(deg)';
                h1.YLabel.String = 'Scheduling Parameter';
                h2.YLabel.String = 'Scheduling Parameter';
                h1.XScale = 'log';
                h2.XScale = 'log';   
            figure;
            [varargout{1:nargout}] = bode(varargin{:});
            if nargout==0, legend(labels{:}), end
        end

        function varargout = bodemag(varargin)
            [varargin,labels] = stdargs(varargin,'freq');
                for i = 1:length(varargin)
                    if (numel(varargin{i}.SamplingGrid) == 1)
                    s =  fieldnames(varargin{i}.SamplingGrid);
                    param_grid{i} = varargin{i}.SamplingGrid.(s{1});
                        for j = 1:size(varargin{i},3)
                            h1 = subplot(1,1,1);
                            [mag,~,w] = bode(varargin{i}(:,:,j,:));
                            hold(h1,'on')
                            h1.ColorOrderIndex = i;
                            plot3(w(:),param_grid{i}(j)*ones(size(w,1),1),mag(:));
                        end
                    end
                end
                h1.View = [45,45];
                h1.XLabel.String = 'Frequency';
                h1.ZLabel.String = 'Magnitude(dB)';
                h1.YLabel.String = 'Scheduling Parameter';
                h1.XScale = 'log';   
            figure;
            [varargout{1:nargout}] = bodemag(varargin{:});
            if nargout==0, legend(labels{:}), end
        end

        function varargout = bodeplot(varargin)
            [varargin,labels] = stdargs(varargin,'freq');
            [varargout{1:nargout}] = bodeplot(varargin{:});
            if nargout==0, legend(labels{:}), end
        end

        function varargout = nyquist(varargin)
            [varargin,labels] = stdargs(varargin,'freq');
            [varargout{1:nargout}] = nyquist(varargin{:});
            if nargout==0, legend(labels{:}), end
        end

        function varargout = nyquistplot(varargin)
            [varargin,labels] = stdargs(varargin,'freq');
            [varargout{1:nargout}] = nyquistplot(varargin{:});
            if nargout==0, legend(labels{:}), end
        end
        
        function varargout = sigma(varargin)
            [varargin,labels] = stdargs(varargin,'freq');
            [varargout{1:nargout}] = sigma(varargin{:});
            if nargout==0, legend(labels{:}), end
        end

        function varargout = freqresp(varargin)
            [varargin,labels] = stdargs(varargin,'freq');
            [varargout{1:nargout}] = freqresp(varargin{:});
            if nargout==0, legend(labels{:}), end
        end

        function varargout = evalfr(varargin)
            varargin = stdargs(varargin,'freq');
            [varargout{1:nargout}] = evalfr(varargin{:});
        end
        
        function varargout = impulse(varargin)
            [varargin,labels] = stdargs(varargin,'time');
            [varargout{1:nargout}] = impulse(varargin{:});
            if nargout==0, legend(labels{:}), end
        end

        function varargout = step(varargin)
            [varargin,labels] = stdargs(varargin,'time');
            [varargout{1:nargout}] = step(varargin{:});
            if nargout==0, legend(labels{:}), end
        end
        
        function varargout = lsim(varargin)
            error('not implemented yet');
        end

        function varargout = pzmap(varargin)
            [varargin,labels] = stdargs(varargin,'linmod');
            [varargout{1:nargout}] = pzmap(varargin{:});
            if nargout==0, legend(labels{:}), end
        end

        function varargout = sim(self,u,varargin)
            models = content(model(self));
            N = length(models);
            t = cell(1,N);
            x = cell(1,N);
            y = cell(1,N);
            labels = cell(1,N);
            for k = 1:length(models)
                [y{k},t{k},x{k}] = sim(ODEmod(models{k}),u,varargin{:});
                labels{k} = models{k}.name;
            end
            
            if nargout == 0
                varargout = {};
                for k = 1:length(self.out)
                    subplot(length(self.out),1,k)
                    yk = cellfun(@(x) x(:,k),y,'un',0);
                    args = [t;yk];
                    plot(args{:});
                end
                subplot(length(self.out),1,1);
                legend(labels{:});
            else
                varargout = {y,t,x};
            end
        end
    end
end

