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

classdef SystemOfSystems < AbstractSystem
    %SYSTEMOFSYSTEMS Summary of this class goes here
    %   Detailed explanation goes here
        
    properties
        connections = {};
    end
    
    methods
        
        [config,K,solver] = solve(config,obj,constr,varargin)
        
        function self = SystemOfSystems(varargin)
            connections_ = {};
            if iscell(varargin{nargin})
                connections_ = varargin{nargin};
                varargin(nargin) = [];
            end
            self = self@AbstractSystem();
            self = add(self,varargin{:});
            self.connections = connections_;
        end
        
        function c = copy(self)
            cont = cellfun(@copy,self.content_,'un',0);
            conn = self.connections;
            c = SystemOfSystems(cont{:},conn);
            c.in = self.in; c.out = self.out;
        end
        
        function self = add(self,varargin)
            assert(all(cellfun(@(x)isa(x,'AbstractSystem'),varargin)),'SystemOfSystems only takes ''AbstractSystems'' as content_.');
            self.content_ = [self.content_,varargin];
            
            in_ = cellfun(@(x)x.in,varargin,'un',0);
            out_ = cellfun(@(x)x.out,varargin,'un',0);
            self.in = [self.in;vertcat(in_{:})];
            self.out = [self.out;vertcat(out_{:})];
            assert(length(unique(self.in))==length(self.in),'Cannot add the same input twice to the same plant');
            assert(length(unique(self.out))==length(self.out),'Cannot add the same output twice to the same plant');
        end
        
        function mod = model(self,varargin)
            mod = unpack(self);
            mod = remove_empty(mod);
            mod = resolve_lincomb(mod);
            mod = resolve_equalio(mod);
            mod = addalias(mod,mod.connections);
            mod = resolve_lft(mod,varargin{:});
        end
        
        function unpacked = unpack(self)
            % UNPACK return an equivalent SoS of SoMs
            
            [systems,conns] = unpack_helper(self);
            unpacked = SystemOfSystems(systems{:},conns);
        end
        
        function [self,restore] = subsystem(self,varargin)
            % SUBSYSTEM return the SoS of which the SoMs retain the
            % specified models
            
            [self,restore] = subsystem_helper(self,varargin{:});
        end
        
        function self = remove_empty(self)
            % REMOVE_EMPTY remove empty systems from the SoS
            
            empty = cellfun(@isempty,self.content_);
            self = SystemOfSystems(self.content_{~empty},self.connections);
        end
    end
    
    methods (Access=protected)
        function disptab(self,tabs)
            disptab@AbstractSystem(self,tabs);
            
            tabs = repmat('\t',[1,numel(tabs)]);
            fprintf([tabs,'\tContent (%d):\n'],length(self.content_))
            if ~isempty(self)
                for k=1:length(self.content_)
                    disptab(self.content_{k},tabs+1);
                end
                fprintf('\n');
            end
        end   
        
        function submod = subsref_signal(self,in,out)
            % SUBSREF_SIGNAL Get part of the system based on signals

            mod = model(self);
            submod = mod(out,in);
        end
    end
        
    methods (Access=private)
        function [systems,connections] = unpack_helper(self)
            [systems,connections] = cellfun(@unpack_helper,self.content_,'un',0);
            systems = horzcat(systems{:});
            connections = vertcat(self.connections,connections{:});
        end
        
        function [self,restore] = subsystem_helper(self,varargin)
            [self.content_,restore] = cellfun(@(x)subsystem_helper(x,varargin{:}),self.content_,'un',0);
            restore = vertcat(restore{:});
        end
                   
        function emptysys = get_empty(self)
            % GET_EMPTY returns all empty systems in the SoS
            
            empty = cellfun(@isempty,self.content_);
            emptysys = self.content_(empty);
        end
                
        function self = resolve_lincomb(self)
            % RESOLVE_LINCOMB add linear combinations of signals to the SoS
            
            % Check for unique signals first: avoid having the same adder
            % twice
            unique_conn = unique(vertcat(self.connections{:}));
            
            for k = 1:numel(unique_conn)
                if islincomb(unique_conn(k))
                    [mod,new_connection] = system(unique_conn(k));
                    self = add(self,mod);
                    self.connections = [self.connections;new_connection];
                end
            end
        end
        
        function self = resolve_equalio(self)
            % RESOLVE_EQUALIO connect equivalent inputs and outputs
            connections_ = Signal.group(self.connections);
            
            for k = 1:length(connections_)
                [membin,~] = ismember(self.in,connections_{k});
                [membout,~] = ismember(self.out,connections_{k});
                membin = find(membin);
                membout = find(membout);
                if (length(membin) == 1) && (length(membout) == 1)
                    self.fbsigin = [self.fbsigin;self.out(membout)];
                    self.fbsigout = [self.fbsigout;self.in(membin)];
                    self.fb = blkdiag(self.fb,1);
                elseif (length(membin) > 1) && (length(membout) <= 1)
                    self.fbsigin = [self.fbsigin;connections_{k}];
                    self.fbsigout = [self.fbsigout;self.in(membin)];
                    self.fb = blkdiag(self.fb,ones(length(membin),1));
                elseif (length(membin) <= 1) && (length(membout) > 1)
                    self.fbsigin = [self.fbsigin;self.out(membout)];
                    self.fbsigout = [self.fbsigout;connections_{k}];
                    self.fb = blkdiag(self.fb,ones(1,length(membout)));
                elseif isempty(membin) || isempty(membout)
                    % ignore this -> pure alias
                else
                    error('InputOutput cluster discovered. Inform the developers.');
                end
            end
        end
        
        function LFT = resolve_lft(self,varargin)
            % RESOLVE_LFT Compute lft of plant and feedback matrix
            
            % 1. do 2-step connection
            [conin,idx_w1,confbout,~] = SystemOfSystems.twostepconn(self.in,self.fbsigout);
            [conout,idx_z1,confbin,idx_w2] = SystemOfSystems.twostepconn(self.out,self.fbsigin);
            idx_z2 = 1:length(self.fbsigout);
            ny = length(conout);
            nu = length(conin);
            
            % 2. compose input/output list for lft
            idx_in = [idx_w1,conin];
            idx_out = [idx_z1,conout];
            idx_fbin = [confbin,idx_w2];
            idx_fbout = [confbout,idx_z2];
            
            % 3. compose content_
            models = cellfun(@(x)model(x,varargin{:}),self.content_,'UniformOutput',false);
            content__ = cellfun(@(x){x},models{1}.content_,'UniformOutput',false);
            for k = 2:length(models)
                content__ = repmat(content__,[length(models{k}.content_),1]);
                for i = 1:size(content__,1)
                    for j = 1:size(content__,2)
                        content__{i,j} = [content__{i,j},models{k}.content_(i)];
                    end
                end
                content__ = content__(:)';
            end
            
            % 4. Safe blockdiagonal
            function [blk,remove] = safeblkdiag(varargin)
                 try
                    blk = blkdiag(varargin{:});
                    names = cellfun(@(x) x.name,varargin,'un',0);
                    blk.name = strjoin(names,'-');
                    remove = false;
                 catch
                     warning('Some models cannot be combined');
                     blk = [];
                     remove = true;
                 end
            end
            [cont,remove] = cellfun(@(x)safeblkdiag(x{:}),content__,'UniformOutput',false);
            cont(cell2mat(remove)) = [];
            
            % 5. compute lft
            content__ = cellfun(@(x)lft(x(idx_out,idx_in),self.fb(idx_fbout,idx_fbin),nu,ny),cont,'UniformOutput',false);
            in_ = [self.in(idx_w1);self.fbsigin(idx_w2)];
            out_ = [self.out(idx_z1);self.fbsigout(idx_z2)];
            
            % 6. add input signals to output
            content__ = cellfun(@(x)lft(x,repmat(eye(length(in_)),[2,1]),length(in_),0),content__,'un',0);
            out_ = [out_;in_];
            LFT = IOSystem(in_,out_);
            LFT = add(LFT,content__{:});
        end
    end
    
    methods (Static,Access=private)
        function [conA,disconA,conB,disconB] = twostepconn(A,B)
            % 0. make variables
            disconA = 1:length(A); disconB = 1:length(B);
            
            % 1. find exact matches
            [~,conexactA,conexactB] = intersect(A,B,true);
            disconA(conexactA) = []; disconB(conexactB) = [];
            A(conexactA) = []; B(conexactB) = [];
            
            % 2. find inexact matches
            [~,conaliasA,conaliasB] = intersect(A,B);
            
            % 3. connections
            conA = [conexactA,disconA(conaliasA)]; conB = [conexactB,disconB(conaliasB)]; % remapped
            disconA(conaliasA) = []; disconB(conaliasB) = [];
        end
    end
    
end

