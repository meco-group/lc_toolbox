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
% terms of @f$\mathcal{H}_\infty@f$ and @f$\mathcal{H}_2@f$norms. 
    
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
        %  varargin : options that are passed to the \c Norm construction,
        %  first part is a \c double defining which norm you're using (\c 2 or
        %  \c Inf)
        %
        % Return values:
        %  n : the norm @type Norm 
            n = Norm(self);
            if nargin > 1
                n = setnorm(n,varargin{1});
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
        % Returns the product of the norms of both input arguments.
        % Overloads MATLAB's \c * operator for Channel objects.
        % 
        % Only of of both input arguments is allowed to be a Channel
        % objects. The product of two channels has no meaning.
        %
        % Parameters: 
        %  tf1 : either a channel, a model (\c numlti or Model) or a
        %  scalar (\c double) 
        %  tf2 : either a channel, a model (\c numlti or Model) or a
        %  scalar (\c double) 
        %
        % Return values:
        %  n : norm of the product of both input arguments @type Norm
            if isnumeric(tf1) && isscalar(tf1)
                n = tf1*Norm(tf2);
            elseif isnumeric(tf2) && isscalar(tf2)
                n = tf2*Norm(tf1);
            else
                n = Norm(tf1,tf2);
            end
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
        % Return values:
        %  sorted : same channel with sorted subsignals @type Channel
            v = arrayfun(@(x)arrayfun(@(y)sum([y.in.UUID y.out.UUID]),x),self);
            [~,i] = sort(v);
            sorted = self(i);
        end

        function [tf] = vertcat(varargin)
        % Concatenates channels vertically. Only possible if the
        % channels have the same input signal(s), e.g. @f$ \left[w
        % \rightarrow z_1 ; w \rightarrow z_2\right]@f$ becomes @f$ w \rightarrow \left[\begin{array}{l} z_1 \\ z_2 \end{array}\right]@f$. 
        % Overloads MATLAB's \c vertcat() for Channel objects.
        % 
        % Parameters:
        %  varargin : list of channels
        %
        % Return values:
        %  tf : channel stacking all outputs in the same order as provided by the user, the input signal remains the same @type Channel
        
            assert(all(cellfun(@(x,y) isequal(x.in,y.in), varargin, circshift(varargin,1))),'Vertical concatenation of channels is only possible if the input signals are equal. Use blkdiag instead.');
            out_ = [];
            for k = 1:length(varargin)
                out_ = [out_;varargin{k}.out];
            end
            if length(varargin)>1
                tf = Channel(varargin{1}.in,out_,['[' strjoin(cellfun(@(x) {x.name}, varargin), '; ') ']']);
            else
                tf = varargin{1};
            end
        end
        
        function [tf] = horzcat(varargin)
        % Concatenates channels horizontally. Only possible if the
        % channels have the same output signal(s), e.g. @f$ \left[w_1
        % \rightarrow z , w_2 \rightarrow z\right]@f$ becomes @f$\left[\begin{array}{l} w_1 \\ w_2 \end{array}\right] \rightarrow z @f$.  
        % Overloads MATLAB's \c horzcat() for Channel objects.
        % 
        % Parameters:
        %  varargin : list of channels
        %
        % Return values:
        %  tf : channel stacking all inputs in the same order as provided by the user, the output signal remains the same @type Channel
        
            assert(all(cellfun(@(x,y) isequal(x.out,y.out), varargin, circshift(varargin,1))),'Horizontal concatenation of channels is only possible if the output signals are equal. Use blkdiag instead.');
            in_ = []; 
            for k = 1:length(varargin)
                in_ = [in_;varargin{k}.in];
            end
            if length(varargin)>1
                tf = Channel(in_,varargin{1}.out,['[' strjoin(cellfun(@(x) {x.name}, varargin), ', ') ']']);
            else
                tf = varargin{1};
            end
        end
        
        function [tf] = blkdiag(varargin)
        % Concatenates channels diagonally, i.e. stacking both input and
        % output signals, e.g. @f$ \text{blkdiag}(w_1
        % \rightarrow z_1, w_2 \rightarrow z_2)@f$ becomes @f$ \left[\begin{array}{l} w_1 \\ w_2 \end{array}\right] \rightarrow \left[\begin{array}{l} z_1 \\ z_2 \end{array}\right]@f$. 
        % Overloads MATLAB's \c blkdiag() for Channel objects.
        % 
        % Parameters:
        %  varargin : list of channels
        %
        % Return values:
        %  tf : channel stacking both all inputs and all outputs @type Channel
        
            in_ = [];
            out_ = [];
            for k = 1:length(varargin)
                in_ = [in_;varargin{k}.in];
                out_ = [out_;varargin{k}.out];
            end
            if length(varargin)>1
                tf = Channel(in_,out_,['blkdiag(' strjoin(cellfun(@(x) {x.name}, varargin), ', ') ')']);
            else
                tf = varargin{1};
            end
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
        
        function [lia,locb] = ismember(self,other)
        % Checks whether the channels of a Channel array are in another
        % Channel array.
        % Overloads MATLAB's \c ismember() for Channel objects.
        % 
        % Parameters:
        %  self: array of channels that is searched for in other @type Channel
        %  other: array of channels that is checked for the presence of channels of self @type Channel
        %
        % Return values:
        %  lia: 1 where the channels of self are found in other @type logical
        %  locb: lowest index values such that self = other(locb) @type logical
        
            assert(isa(self,'Channel') && isa(other,'Channel'), 'ismember only works on two Channel objects.');
            if size(self,1)==1 && size(other,1)==1
            elseif size(self,1)==1 && size(other,2)==1
                other = other';
            elseif size(self,2)==1 && size(other,1)==1
                self = self';
            elseif size(self,2)==1 && size(other,2)==1
                self = self'; other = other';
            else
                error('ismember only works on vectors of Channel objects.');
            end

            lia = zeros(size(self));
            locb = zeros(size(self));
            for i=1:length(self)
                found = false; j = 1;
                while j <= length(other) && ~found
                    if self(i)==other(j)
                        found = true;
                        lia(i) = 1;
                        locb(i) = j;
                    end
                    j=j+1;
                end
            end
            
        end
        
    end
    
end

