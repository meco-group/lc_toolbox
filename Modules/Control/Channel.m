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

classdef(InferiorClasses = {?DysSys,?zpk,?ss,?tf}) Channel
% Channel objects are used to represent an input-output relation between
% signals. These could then be used to specify performance objectives in
% terms of H-infinity and H-2 norms. 
    
    properties
        in;         % the input signal of the channel @type Signal
        out;        % the output signal of the channel @type Signal
        name = '';  % the name of the channel @type char
    end
    
    methods
        function[obj] = Channel(varargin)
        % Constructor for Channel objects.
        % 
        % Parameters:
        %  varargin : should have on of these three structures:
        %   -# an existing Channel, followed by a name (\c char) you want to
        %   assign to it
        %   -# an input Signal, followed by an output Signal
        %   -# same as option 2, followed by a name (\c char) as last input argument
        %
        % Return values:
        %  self: the channel @type Channel 
            if(nargin>=2)
                if(isa(varargin{1},'Channel'))
                    obj = varargin{1};
                    obj.name = varargin{2};
                else
                    obj.in = varargin{1};
                    obj.out = varargin{2};
                end
            end
            if(nargin>=3)
                obj.name = varargin{3};
            end
        end
        
        function n = norm(self,varargin)
        % Defines the norm of the channel.
        % 
        % Parameters:
        %  varargin : options that are passed to the \c norm construction,
        %  first part is a \c double defining which norm you're using (\c 2 or
        %  \c Inf)
        %
        % Return values:
        %  n : the norm @type Norm 
            n = Norm(self);
            if nargin > 1
                n = norm(n,varargin{1});
            end
        end
                
        function s = csize(self)
        % Returns the size of the channel. 
        %
        % Return values:
        %  s : [number of outputs, number of inputs] @type double 
            s = [length(self.out),length(self.in)];
        end
        
        function b = issiso(self)
        % Checks whether the channel is SISO.
        %
        % Return values:
        %  b : boolean reflecting whether the channel has only 1 input signal and 1 output signal @type logical 
            b = arrayfun(@(x)all(csize(x)==1),self);
        end
        
        function n = mtimes(tf1,tf2)
        % Return the product of the norms of both input arguments.
        % Overloads MATLAB's \c * operator for Channel objects.
        % 
        % Only of of both input arguments is allowed to be a Channel
        % objects. The product of two channels has no meaning.
        %
        % Parameters: 
        %  tf1 : either a channel, a standard MATLAB model (\c numlti) or a
        %  scalar (\c double) 
        %  tf2 : either a channel, a standard MATLAB model (\c numlti) or a
        %  scalar (\c double) 
        %
        % Return values:
        %  n : norm of the product of both input arguments @type Norm
            n = Norm(tf1,tf2);
        end
        
        function b = eq(tf1,tf2)
        % Checks whether two channels are the same. 
        % Overloads MATLAB's \c == operator for Channel objects.
        %
        % @note \c eq uses Signal::isequal to check whether or not the
        % input and output channels are the same and thus does not consider
        % signals that might be aliases (de facto the same). 
        %
        % Parameters: 
        %  tf1 : first channel @type Channel
        %  tf2 : second channel @type Channel
        %
        % Return values:
        %  b : boolean reflecting whether both channels are equal
        %  @type logical
            assert(numel(tf1) == numel(tf2),'eq only accepts equally sized arrays');
            if numel(tf1) == 1
                b = false;
                if all(csize(tf1) == csize(tf2))
                    bin = all(arrayfun(@(x,y)isequal(x,y),tf1.in,tf2.in));
                    bout = all(arrayfun(@(x,y)isequal(x,y),tf1.out,tf2.out));
                    b = bin && bout;
                end
            else
                b = arrayfun(@(x,y)eq(x,y),tf1,tf2);
            end
        end
        
        function b = ne(tf1,tf2)
        % Checks whether two channels are \b not the same. 
        % Overloads MATLAB's \c ~= operator for Channel objects.
        %
        % Parameters: 
        %  tf1 : first channel @type Channel
        %  tf2 : second channel @type Channel
        %
        % Return values:
        %  b : boolean reflecting whether both channels are not equal
        %  @type logical
            b = ~eq(tf1,tf2);
        end
        
        function sorted = sort(self)
        % Sorts the subsignals of the input and output signal based on their identifiers.  
        % Overloads MATLAB's \c sort for Channel objects.
        %
        % Parameters:
        %  tf1 : first channel @type Channel
        %  tf2 : second channel @type Channel
        %
        % Return values:
        %  sorted : same channel with sorted subsignals @type Channel
            v = arrayfun(@(x)arrayfun(@(y)sum([y.in.UUID y.out.UUID]),x),self);
            [~,i] = sort(v);
            sorted = self(i);
        end
        
        function [b,io] = abstractEquality(tf1,tf2)
            if((length(tf1.in)~=length(tf2.in))||(length(tf1.out)~=length(tf2.out)))
                b = false;
                io = [tf1.in;...
                      tf1.out]; 
            else
                s = [tf1.in  tf2.in;...
                     tf1.out tf2.out];            
                b = isequal(s(:,1),s(:,2));
                io = s(:,1);
            end
        end
        
        function [tf] = vertcat(varargin)
        % Concatenates two Channels vertically.
        % Overloads MATLAB's \c vertcat() for Channel objects.
        % 
        % Parameters:
        %  varargin : list of channels
        %
        % Return values:
        %  tf : channel stacking all inputs and outputs of the channels in the same order as provided by the user @type Channel
            in = []; out = [];
            for k = 1:length(varargin)
                in = [in;varargin{k}.in];
                out = [out;varargin{k}.out];
            end
            tf = Channel(in,out,varargin{1}.name);
        end
        
        function varargout = subsref(self,s)
        % Implements MATLAB subscripting for channels.
        % Overloads MATLAB's \c subsref() for Channel objects.
        % 
        % Parameters:
        %  s : \c substruct describing the subscript format
        %
        % Return values:
        %  varargout : depending on the subscript a subchannel (\c Channel)
        %  or a system (\c AbstractSystem)
            if strcmp(s(1).type,'()') && ~(length(s(1).subs)==1 && strcmp(s(1).subs{1},':'))
                switch(length(s(1).subs))
                    case 0
                        varargout = {self};
                    case 1
                        if isa(s(1).subs{1},'AbstractSystem') %model access
                            varargout = {subsref(s(1).subs{1},substruct('()',{self.out,self.in}))};
                        else
                            varargout = {builtin('subsref',self,s(1))};
                        end
                    case 2
                        siz = csize(self);
                        if s(1).subs{1}==':';s(1).subs{1}=1:siz(1);end
                        if s(1).subs{2}==':';s(1).subs{2}=1:siz(2);end
                        assert(all(cellfun(@isnumeric,s(1).subs)),'only numeric index');
                        assert(all(s(1).subs{1}<=siz(1)) && all(s(1).subs{2}<=siz(2)),'Index exceeds system dimensions.');
                        
                        varargout = {Channel(self.in(s(1).subs{2}),self.out(s(1).subs{1}))};
                end
                if length(s) > 1
                    varargout = {builtin('subsref',varargout{1},s(2:end))};
                end
            else
                varargout = {builtin('subsref',self,s(:))};
            end
        end
    end
    
end

