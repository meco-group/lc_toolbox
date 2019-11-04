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

classdef SystemOfModels < AbstractSystem
    %SYSTEMOFMODELS Summary of this class goes here
    %   Detailed explanation goes here
    
    methods
        function self = SystemOfModels(varargin)
            if all(cellfun(@(x)isa(x,'Model'),varargin))
                siz = {size(varargin{1},2),size(varargin{1},1)};
            else
                siz = varargin;
            end
            self = self@AbstractSystem(siz{:});
            if isa(varargin{1},'Model')
                self = add(self,varargin{:});
            end
        end
        
        function c = copy(self)
            c = SystemOfModels(self.in,self.out);
            c.content_ = self.content_;
        end
        
        function self = add(self,varargin)
            assert(all(cellfun(@(x)isa(x,'Model'),varargin)),'SystemOfModels only takes ''Model'' as content_.');
            assert(all(cellfun(@(x)all(size(x)==size(self)),varargin)),'Dimensions of added model mismatch');
            for k = 1:length(varargin)
                if isempty(varargin{k}.name)
                    varargin{k}.name = ['model' num2str(length(self.content_)+k)];
                end
            end
            self.content_ = [self.content_,varargin];
        end
        
        function mod = model(self,varargin)
            mod = copy(self);
            mod = subsystem_helper(mod,varargin{:});
        end
        
        function restore(self,removed,removedindex)
            if nargin == 2
                removedindex = 1:length(removed);
            end
            contentindex = 1:(length(removed)+length(self.content_));
            contentindex(removedindex) = [];
            [~,order] = sort([contentindex,removedindex]);
            
            self.add(removed{:});
            self.content_ = self.content_(order);
        end
    end
    
    % This method should be only visible to the SystemOfSystems
    % class (friend access, matlab styled?)
    methods
        function [systems,connections] = unpack_helper(self)
            systems = {self};
            connections = {};
        end
        
        function [self,restore] = subsystem_helper(self,varargin)
            restore = [];
            if nargin > 1
                if all(cellfun(@ischar,varargin))
                    match = cellfun(@(x)any(strcmp(x.name,varargin)),self.content_);
                    index = 1:numod(self);
                    keep = index(match);
                else
                    assert(isnumeric(varargin{1}),'Input should be char or numeric');
                    keep = varargin{1};
                end

                if ~isempty(keep) || isempty(varargin{1})
                    rmindex = 1:numod(self);
                    rmindex(keep) = [];
                    restore = struct('var',self,'removed',{self.content_(rmindex)},'index',{rmindex});
                    self.content_(rmindex) = [];
                end
            end
        end
    end
    
    methods (Access=protected)
        function disptab(self,tabs)
            disptab@AbstractSystem(self,tabs);
            
            tabs = repmat('\t',[1,numel(tabs)]);
            fprintf([tabs,'\tContent (%d):\n'],length(self.content_))
            if ~isempty(self)
                fprintf([tabs '%s'],self.content_{1}.name);
                for k = 2:length(self.content_)
                    fprintf(', %s',self.content_{k}.name);
                end
                fprintf('\n');
            end
        end
        
        function submod = subsref_signal(self,in,out)
            function str = signalsNotFound(sig)
                sig = listaliases(sig);
                if length(sig)>1
                    l = sprintf('%i, ', sig(2:end)); l = l(1:end-2);
                    str = [num2str(sig(1)) ' (also known as ' l ')']; 
                else
                    str = num2str(sig);
                end
            end
            
            [mtcin,idxin] = ismember(in,self.in);
            [mtcout,idxout] = ismember(out,self.out);
            
            if all(mtcin) && all(mtcout)
                content__ = cellfun(@(x)x(idxout,idxin),self.content_,'UniformOutput',false);
                submod = SystemOfModels(self.in(idxin(idxin~=0)),self.out(idxout(idxout~=0)));
                submod = add(submod,content__{:});
            else
                lcin = in(~mtcin);
                connin = {};
                if ~all(mtcin)
                    for k = 1:length(lcin) %connect inputs to make the linear comb
                        s = lcin(k);
                        lcin(k).multipliers = 1;
                        lcin(k).components = [];
                        if ~all(size(s.components)>=1)
                            error(['The following signal was not found in the system: ' signalsNotFound(s) '. You might have to check your plant, as this usually occurs when accidentally closing a loop: your exogenous signals then become internal signals.']);
                        end
                        if ~all(ismember(s.components,self.in))
                            nip = s.components(~ismember(s.components,self.in)); 
                            for i=1:length(nip)
                                err{i} = signalsNotFound(nip(i));
                            end
                            error(['The following signal(s), that are part of a linear combination signal you requested, were not members of the plant: ' strjoin(err,', ') '.']);
                        end
                        connin_ = arrayfun(@(m,c) eq(c,m*lcin(k)),s.multipliers(:),s.components(:),'un',0);
                        connin = vertcat(connin,connin_{:});
                    end
                    in(~mtcin) = lcin;
                end
            
                lcout = out(~mtcout);
                connout = {};
                if ~all(mtcout)
                    for k = 1:length(lcout) %connect outputs to make linear comb
                        s = lcout(k);
                        lcout(k).multipliers = 1;
                        lcout(k).components = [];
                        if ~all(size(s.components)>=1)
                            error(['The following signal was not found in the system: ' signalsNotFound(s) '. You might have to check your plant, as this usually occurs when accidentally closing a loop: your exogenous signals then become internal signals.']);
                        end
                        if ~all(ismember(s.components,[self.out;self.in]))
                            nip = s.components(~ismember(s.components,[self.out;self.in])); 
                            for i=1:length(nip)
                                err{i} = signalsNotFound(nip(i));
                            end
                            error(['The following signal(s), that are part of a linear combination signal you requested, were not members of the plant: ' strjoin(err,', ') '.']);
                        end
                        connout_ = arrayfun(@(m,c) eq(m*c,lcout(k)),s.multipliers,s.components,'un',0);
                        connout = vertcat(connout,connout_{:});
                    end
                end
            
                conn = [connin;connout];
                P = IOSystem(self,conn);
                P = model(P);
                submod = subsref_signal(P,in,out);
            end
        end
    end
end

