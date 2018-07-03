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

classdef (InferiorClasses = {?zpk,?tf,?ss,?frd}) Model
    %MODEL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        name = '';
    end
    
    methods (Abstract)
        nin(self);
        nout(self);
        std(self);
        blkdiag(varargin);
        lft(self,other,nu,ny);
    end
    
    methods (Abstract, Access=protected)
        submodel(self,idxout,idxin);
    end
    
    methods
        function self = Model()
            self.name = '';
        end
        
        function self = setName(self,name)
            self.name = name;
        end
        
        function [varargout] = size(self,dim)
            siz = [nout(self),nin(self)];

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
            % SUBSREF reimplement subsref for systems subindexing            
            if strcmp(s(1).type,'()')
                siz = size(self);
                assert(~(isa(s(1).subs{1},'Signal') || isa(s(1).subs{2},'Signal')),'Cannot index a model using signals. Signals refer only to systems');
                if s(1).subs{1}==':';s(1).subs{1}=1:siz(1);end
                if s(1).subs{2}==':';s(1).subs{2}=1:siz(2);end
                assert(all(cellfun(@isnumeric,s(1).subs)),'only numeric index');
                assert(all(s(1).subs{1}<=siz(1)) && all(s(1).subs{2}<=siz(2)),'Index exceeds system dimensions.');
                varargout = {submodel(self,s(1).subs{1},s(1).subs{2})};
                
                if length(s)>1
                    varargout = {builtin('subsref',varargout{1},s(2:end))};
                end
            else
                varargout = {builtin('subsref',self,s(:))};
            end
        end
        
        function sum = plus(self,other)
            assert(isa(self,'Model') || isa(self,'numlti') || isnumeric(self),'Both terms should be models.');
            assert(isa(other,'Model') || isa(other,'numlti') || isnumeric(other),'Both terms should be models.');
            assert(all(size(self) == size(other)),'Model I/O dimensions must agree.');
            
            one = IOSystem(fromstd(self));
            two = IOSystem(fromstd(other));
            output = Signal(size(one,1));
            connections = [one.in == two.in ; output == one.out+two.out];
            
            sum = IOSystem(one,two,connections);
            sum = sum.model();
            sum = sum(output,one.in);
            sum = sum.content(1);
        end
        
        function difference = minus(self,other)
            difference = plus(self,uminus(other));
        end
        
        function opposite = uminus(self)
        	opposite = SSmod(-eye(size(self,1)))*self;
        end
        
        function product = mtimes(self,other)
            
            if isa(self,'Channel') || isa(other,'Channel')
                product = Norm(self,other);
            else
            
            assert(isa(self,'Model') || isa(self,'numlti') || isnumeric(self),'Both terms should be models.');
            assert(isa(other,'Model') || isa(other,'numlti') || isnumeric(other),'Both terms should be models.');
            
            if isnumeric(self) && isscalar(self), self = self*eye(size(other,1)); end
            if isnumeric(other) && isscalar(other), other = other*eye(size(self,2)); end
            assert(size(self,2) == size(other,1),'Number of outputs of the first model does not equal number of inputs of the second model.');

            one = IOSystem(fromstd(self));
            two = IOSystem(fromstd(other));
            connections = [two.out == one.in];

            product = IOSystem(one,two,connections);
            product = product.model();
            product = product(one.out,two.in);
            product = product.content(1);
            end
        end
       
       
        function cat = vertcat(self,varargin)
            
            assert(nargin >= 2,'vertcat needs at least two arguments.'); 
            list = {self, varargin{:}};
            assert(range(cellfun(@(x) size(x,2),list))==0,'Number of inputs should be equal.');
            assert(all(cellfun(@(x) (isa(x,'Model') || isa(x,'numlti') || isnumeric(x)),list)),'You can only concatenate models or models and constant matrices.');

            connections = cell(length(list)-1,1);
            if isa(list{1},'numlti'); list{1} = fromstd(list{1}); end;
            if isnumeric(list{1}); list{1} = SSmod(list{1}); end;
            list{1} = IOSystem(list{1});
            for i = 2:length(list)
                if isa(list{i},'numlti'); list{i} = fromstd(list{i}); end;
                if isnumeric(list{i}); list{i} = SSmod(list{i}); end;
                list{i} = IOSystem(list{i});
                connections{i-1} = list{i-1}.in == list{i}.in;
            end
            
            listout = cellfun(@(x) x.out,list,'UniformOutput',false);
            cat = IOSystem(list{:},vertcat(connections{:}));
            cat = cat.model();
            cat = cat(vertcat(listout{:}),list{1}.in);
            cat = cat.content(1);
                    
        end
        
        function cat = horzcat(self,varargin)
            
            assert(nargin >= 2,'horzcat needs at least two arguments.'); 
            list = {self, varargin{:}};
            assert(range(cellfun(@(x) size(x,1),list))==0,'Number of outputs should be equal.');
            assert(all(cellfun(@(x) (isa(x,'Model') || isa(x,'numlti') || isnumeric(x)),list)),'You can only concatenate models or models and constant matrices.');

            connections = cell(length(list)-1,1);
            if isa(list{1},'numlti'); list{1} = fromstd(list{1}); end;
            if isnumeric(list{1}); list{1} = SSmod(list{1}); end;
            list{1} = IOSystem(list{1});
            for i = 2:length(list)
                if isa(list{i},'numlti'); list{i} = fromstd(list{i}); end;
                if isnumeric(list{i}); list{i} = SSmod(list{i}); end;
                list{i} = IOSystem(list{i});
                connections{i-1} = list{i-1}.out == list{i}.out;
            end
            
            listin = cellfun(@(x) x.in,list,'UniformOutput',false);
            cat = IOSystem(list{:},vertcat(connections{:}));
            cat = cat.model();
            cat = cat(list{1}.out,vertcat(listin{:}));
            cat = cat.content(1);

        end
        
    end
    
    % Plotting
    methods (Access=public)
        function varargout = bode(varargin)
            varargin = stdargs(varargin,'freq');
            [varargout{1:nargout}] = bode(varargin{:});
        end
        
        function varargout = bodemag(varargin)
            varargin = stdargs(varargin,'freq');
            [varargout{1:nargout}] = bodemag(varargin{:});
        end
        
        function varargout = bodeplot(varargin)
            varargin = stdargs(varargin,'freq');
            [varargout{1:nargout}] = bodeplot(varargin{:});
        end
        
        function varargout = nyquist(varargin)
            varargin = stdargs(varargin,'freq');
            [varargout{1:nargout}] = nyquist(varargin{:});
        end
        
        function varargout = nyquistplot(varargin)
            varargin = stdargs(varargin,'freq');
            [varargout{1:nargout}] = nyquistplot(varargin{:});
        end
        
        function varargout = sigma(varargin)
            varargin = stdargs(varargin,'freq');
            [varargout{1:nargout}] = sigma(varargin{:});
        end
        
        function varargout = freqresp(varargin)
            varargin = stdargs(varargin,'freq');
            [varargout{1:nargout}] = freqresp(varargin{:});
        end
        
        function varargout = evalfr(varargin)
            varargin = stdargs(varargin,'freq');
            [varargout{1:nargout}] = evalfr(varargin{:});
        end
        
        function varargout = impulse(varargin)
            varargin = stdargs(varargin,'time');
            [varargout{1:nargout}] = impulse(varargin{:});
        end
                
        function varargout = step(varargin)
            varargin = stdargs(varargin,'time');
            [varargout{1:nargout}] = step(varargin{:});
        end
             
        function varargout = pzmap(varargin)
            varargin = stdargs(varargin,'linmod');
            [varargout{1:nargout}] = pzmap(varargin{:});
        end   
    end
end

