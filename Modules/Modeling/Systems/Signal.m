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

classdef(InferiorClasses = {?zpk,?ss,?tf}) Signal
% Signal objects are used to represent signals and allow to define
% intuitive aliases, input-output relations (channels) and connections. 
%
% Signals are only used as vectors within LCToolbox. Most methods in this
% class are therefore only applicable for ... x 1 Signal objects. 
    
    properties
        UUID;                       % a unique identifier for the signal in the MATLAB workspace @type double
        components = zeros(1,0);    % components of the linear combination if the signal is a linear combination of multiple other signals @type Signal
        multipliers = 1;            % scaling factor of every signal if a signal is a linear combination of multiple other signals @type double
        alias = [];                 % aliases of the signal that are known @type double
    end
    
    methods (Static = true)
        function out = newUUID()
        % Creates a new unique identifier. 
        %
        % Return values:
        %  out : a unique identifier @type double
            persistent UUID;
            if(~isempty(UUID))
                UUID = UUID+1;
            else
                UUID = 1;
            end
            out = UUID;
        end
    end
    
    methods (Static)
        function self = LinearCombination(signals,multipliers)
        % Creates a new signal from a linear combination of multiple other signals.
        %
        % If a subsignal is already a linear combination on itself, the new
        % linear combination will also unpack this subsignal into its
        % components. 
        % 
        % Parameters:
        %  signals : an 1 x N Signal object, N representing the number of subsignals to be combined @type Signal
        %  multipliers : vector of size N representing the weighting factors to be included in the linear combination @type double
        %
        % Return values:
        %  self : the linear combination (1 x 1) @type Signal
        
            [N,M] = size(signals);
            self = Signal(N);
            
            % Construct overall list of components and multipliers
            components_ = cell(1,M);
            multipliers_ = cell(1,M);
            for k = 1:M
                components_{k} = vertcat(signals(:,k).components);
                if isempty(components_{k})
                    components_{k} = signals(:,k);
                end
                multipliers_{k} = vertcat(signals(:,k).multipliers)*multipliers(k);
            end
            components_ = horzcat(components_{:});
            multipliers_ = horzcat(multipliers_{:});
            
            % Construct new signal for linear combination
            for k = 1:N
                self(k).components = components_(k,:);
                self(k).multipliers = multipliers_(k,:);
            end
        end
        
        function groups = group(connections)
        % Implements connections by assigning aliases to signals.
        % 
        % Parameters:
        %  connections : two-column cell representing the connections; every 
        %  signal of the first column equals the signal in the second 
        %  column @type cell
        %
        % Return values:
        %  group : cell with signals and their aliases that implicitly
        %  define the connections @type cell
            if ~isempty(connections)
                groups = {connections{1,1}};
                for k = 2:size(connections,1)
                    i = 1;
                    while i <= length(groups)
                        if isalias(connections{k,1},groups{i})
                            groups{i} = addalias(groups{i},connections{k,1});
                            break;
                        end
                        i = i+1;
                    end
                    if i > length(groups)
                        groups{i} = connections{k,1};
                    end
                end
            else
                groups = {};
            end
        end
    end
    
    methods          
        function self = Signal(varargin)
        % Constructor for Signal objects.
        % 
        % Parameters:
        %  varargin : may contain the dimension of the signal (\c double); 1 by default. 
        %
        % Return values:
        %  self: the signal @type Signal 
        
            switch nargin
                case 0
                    self.UUID = Signal.newUUID();
                case 1
                    self(varargin{1},1) = Signal();
                    for k = 1:varargin{1} 
                        self(k,1) = Signal();
                    end
            end
        end
                
        %% Implemented relational operators on signals 
        
        function assertsize(s1,s2)
        % Checks whether two signals have the same dimension.
        % 
        % Parameters:
        %  s1 : first signal @type Signal
        %  s2 : second signal @type Signal
            assert(length(s1) == length(s2),'cannot operate on signals of different length')
        end
        
        function connection = connect(s1,s2)
        % Connects two signals by defining each other as aliases.
        % 
        % Parameters:
        %  s1 : first signal @type Signal
        %  s2 : second signal @type Signal
        %
        % Return values:
        %  connection : two-column cell representing the connections; the 
        %  signal of the first column equals the signal in the second 
        %  column  @type cell
            assertsize(s1,s2);
            [s1,s2] = arrayfun(@(x,y)addalias(x,y),s1,s2,'UniformOutput',false);
            connection = num2cell([vertcat(s1{:}),vertcat(s2{:})]);
        end
        
        function s = plus(s1,s2)
        % Adds two signals to each other.
        % Overloads MATLAB's \c + operator for Signal objects.
        % 
        % Parameters:
        %  s1 : first signal @type Signal
        %  s2 : second signal @type Signal
        %
        % Return values:
        %  s : sum of \c s1 and \c s2 @type Signal       
            s = Signal.LinearCombination([s1,s2],[1,1]);
        end
        
        function s = minus(s1,s2)
        % Substracts a signal from another one.
        % Overloads MATLAB's \c - operator for Signal objects.
        % 
        % Parameters:
        %  s1 : first signal @type Signal
        %  s2 : second signal @type Signal
        %
        % Return values:
        %  s : difference between \c s1 and \c s2 @type Signal  
        s = Signal.LinearCombination([s1,s2],[1,-1]);
        end
        
        function s = uminus(s1)
        % Inverts the sign of a signal.
        % Overloads MATLAB's \c - operator for Signal objects.
        % 
        % Parameters:
        %  s1 : signal @type Signal
        %
        % Return values:
        %  s : signal \c s1 with the opposite sign @type Signal
        s = Signal.LinearCombination([s1],[-1]);
        end
        
        function s = mtimes(a,b)
        % Allows to multiply a constant (scalar) with a Signal object.
        % Overloads MATLAB's \c * operator for Signal objects.
        % 
        % Only one of the input arguments is allowed to be a Signal object.
        %
        % Parameters:
        %  a : signal or numeric value
        %  b : signal or numeric value
        %
        % Return values:
        %  s : weighted signal @type Signal
        
            if isnumeric(a)
                m = a;
                s1 = b;
            else
                m = b;
                s1 = a;
            end
            s = Signal.LinearCombination([s1],[m]);
        end
        
        function b = islincomb(self)
        % Checks whether the signal is a linear combination of other
        % signals.
        %
        % Return values:
        %  b : boolean reflecting a signal that is a linear combination
        %  @type logical
            b = arrayfun(@(x)~isempty(x.components),self);
        end     
        
        function [x,y] = addalias(x,y)
        % Adds aliases to two signals that are equal. 
        %
        % Parameters:
        %  x : first signal @type Signal
        %  y : second signal @type Signal
        %
        % Return values:
        %  x : first signal with second signal appended as alias @type Signal
        %  y : second signal with first signal appended as alias @type Signal
            x.alias = unique([x.alias,y,y.alias]);
            if all(size(y)==1)
                y.alias = unique([y.alias,x,x.alias]);
            end
        end
        
        %% Relational operators
        function conn = eq(s1,s2)
        % States that two signals are equal and therefore connects them.
        % Overloads MATLAB's \c == operator for Signal objects.
        %
        % Parameters:
        %  s1 : first signal @type Signal
        %  s2 : second signal @type Signal
        %
        % Return values:
        %  conn : two-column cell representing the connections; the 
        %  signal of the first column equals the signal in the second 
        %  column  @type cell
            conn = connect(s1,s2);
        end
        
        function b = ne(s1,s2)
        % Checks whether two signals are \b not equal.
        % Overloads MATLAB's \c ~= operator for Signal objects.
        %
        % Parameters:
        %  s1 : first signal @type Signal
        %  s2 : second signal @type Signal
        %
        % Return values:
        %  b : boolean reflecting whether or not the signals are the same
        %  @type logical
            b = ~isequal(s1,s2);
        end
        
        function b = lt(s1,s2)
        % Checks whether the identifier of the first signal is smaller than
        % the identifier of the second signal.
        % Overloads MATLAB's \c < operator for Signal objects.
        %
        % Parameters:
        %  s1 : first signal @type Signal
        %  s2 : second signal @type Signal
        %
        % Return values:
        %  b : boolean reflecting whether or not the identifier of the first 
        % signal is smaller than the identifier of the second signal @type logical
            b = lt(s1.UUID,s2.UUID);
        end
        
        function b = gt(s1,s2)
        % Checks whether the identifier of the first signal is larger than
        % the identifier of the second signal.
        % Overloads MATLAB's \c > operator for Signal objects.
        %
        % Parameters:
        %  s1 : first signal @type Signal
        %  s2 : second signal @type Signal
        %
        % Return values:
        %  b : boolean reflecting whether or not the identifier of the first 
        % signal is larger than the identifier of the second signal @type logical
            b = gt(s1.UUID,s2.UUID);
        end
        
        function b = le(s1,s2)
        % Checks whether the identifier of the first signal is smaller than
        % or equal to the identifier of the second signal.
        % Overloads MATLAB's \c <= operator for Signal objects.
        %
        % Parameters:
        %  s1 : first signal @type Signal
        %  s2 : second signal @type Signal
        %
        % Return values:
        %  b : boolean reflecting whether or not the identifier of the first 
        % signal is smaller than or equal to the identifier of the second signal @type logical
            b = ~gt(s1,s2);
        end
        
        function b = ge(s1,s2)
        % Checks whether the identifier of the first signal is larger than
        % or equal to the identifier of the second signal.
        % Overloads MATLAB's \c >= operator for Signal objects.
        %
        % Parameters:
        %  s1 : first signal @type Signal
        %  s2 : second signal @type Signal
        %
        % Return values:
        %  b : boolean reflecting whether or not the identifier of the first 
        % signal is larger than or equal to the identifier of the second signal @type logical
            b = ~lt(s1,s2);
        end
        
        function [Y,I] = sort(X, varargin)
        % Sorts the subsignals based on their identifiers.
        % Overloads MATLAB's \c sort() for Signal objects.
        %
        % Parameters:
        %  varargin : sorting settings, see MATLAB's \c sort
        %
        % Return values:
        %  Y : signal with sorted subsignals @type Signal
        %  I : indices for original signal to obtain sorted list \c Y @type double
        
            x = reshape([X.UUID],size(X));
            [~,I] = sort(x,varargin{:});
            Y = X(I);
        end
        
        function b = isequal(s1,s2)
        % Checks whether two signals are equal.
        % Overloads MATLAB's \c isequal() for Signal objects.
        %
        % @note \c isequal checks for the same size and the same
        % identifiers. Signals that are aliases, and thus de facto equal, 
        % are \b not considered as equal. 
        %
        % Parameters:
        %  s1 : first signal @type Signal
        %  s2 : second signal @type Signal
        %
        % Return values:
        %  b : boolean reflecting whether or not both signals are the same
        %  @type logical
            b = false;
            if isa(s1,'Signal') && isa(s2,'Signal') && (length(s1)==length(s2))
                b = all([s1.UUID] == [s2.UUID]);
            end
        end
        
        function b = isalias(s1,s2)
        % Checks whether the signal is an alias of the other signal. 
        %
        % @note \c isalias is not to be confused with \c isequal. 
        %
        % Parameters:
        %  s1 : first signal @type Signal
        %  s2 : second signal @type Signal
        %
        % Return values:
        %  b : boolean reflecting whether or not the signals are aliases
        %  @type logical
            if length(s1) == length(s2)
                if length(s1) == 1
                    z1 = [s1,s1.alias];
                    z2 = [s2,s2.alias];
                    b = any(arrayfun(@(y)any(arrayfun(@(x)isequal(x,y),z1)),z2));
                else
                    b = all(arrayfun(@(x,y)isalias(x,y),s1,s2));                
                end
            else
                b = false;
            end
        end
        
        function [LIA,LOCB] = ismember(a,b,exact)
        % Checks whether the subsignal of one Signal object are present in
        % the subsignals of another Signal object. 
        % Overloads MATLAB's \c ismember() for Signal objects.
        %
        % Parameters:
        %  a : first signal @type Signal
        %  b : second signal @type Signal
        %  exact : if true, aliases are not considered as equal signals @type logical
        %
        % Return values:
        %  LIA : boolean reflecting whether the subsignals of \c a were found
        %  in \c b @type logical
        %  LOCB : indices of \c b where the values of \a can be found @type
        %  double
            if nargin < 3, exact = false; end
            if exact, compare = @isequal;
            else compare = @isalias; end
            
            LIA = zeros(1,length(a));
            LOCB = zeros(1,length(a));
            for ia = 1:length(a)
                for ib = 1:length(b)
                    if compare(a(ia),b(ib))
                        LIA(ia) = true;
                        LOCB(ia) = ib;
                        break;
                    end
                end
            end
            LIA = logical(LIA);
        end
        
        function [C,IA,IC] = unique(A,varargin)
        % Removes duplicates (not including aliases!) that are present in
        % subsignals.
        % Overloads MATLAB's \c unique() for Signal objects.
        %
        % Parameters:
        %  A : the signal @type Signal
        %  varargin: sorting settings (\c unique)
        %
        % Return values:
        %  C : signal containing the same subsignals as the input but
        %  sorted and without repetitions
        %  IA : indices such that @code C = A(IA) @endcode @type double
        %  IC : indices such that @code A = C(IC) @endcode @type double
%             keyboard
            [~,IA,IC] = unique([A.UUID],varargin{:});
            C = A(IA);
        end
        
        function [C,IA,IB] = intersect(A,B,exact)
        % Returns the set of subsignals that are part of both signals. 
        % Overloads MATLAB's \c intersect() for Signal objects.
        %
        % Parameters:
        %  A : first signal @type Signal
        %  B : second signal @type Signal
        %  exact : if true, aliases are not considered as equal signals @type logical
        %
        % Return values:
        %  C : signal containing the subsignals that are part of both
        %  signals @type Signal
        %  IA : indices such that @code C = A(IA) @endcode @type double
        %  IB : indices such that @code C = B(IB) @endcode @type double
            if nargin < 3, exact = false; end
            
            A = transpose(A(:));
            B = transpose(B(:));
            [U,LIA] = unique(A);
            [LIU,IB] = ismember(U,B,exact);
            C = U(LIU);
            IA = LIA(LIU);
            IA = transpose(IA(:));
            IB(IB==0) = [];
        end
        
        function [C,IA] = setdiff(A,B,exact)
        % Returns the set of subsignals of the first
        % signal that are not in the second signal. 
        % Overloads MATLAB's \c setdiff() for Signal objects.
        %
        % Parameters:
        %  A : first signal @type Signal
        %  B : second signal @type Signal
        %  exact : if true, aliases are not considered as equal signals @type logical
        %
        % Return values:
        %  C : signal containing the subsignals that are part of both
        %  signals @type Signal
        %  IA : indices such that @code C = A(IA) @endcode @type double
            if nargin < 3, exact = false; end
            
            A = transpose(A(:));
            B = transpose(B(:));
            [LIA,~] = ismember(A,B,exact);
            IA = 1:length(A);
            IA(LIA) = [];
            C = A(IA);
        end
        
        %% Implement other stuff..
        
        function TF = mrdivide(s1,s2)
        % Creates a channel from the first signal to the second
        % signal.
        % Overloads MATLAB's \c '/' operator for Signal objects.
        %
        % Parameters:
        %  s1 : first signal @type Signal
        %  s2 : second signal @type Signal
        %
        % Return values:
        %  TF : channel with \c s2 as input and \c s1 as output
        TF = Channel(s2,s1);
        end
                
        function [s] = selection(x,m)
        % Returns a selection matrix to select specific subsignals. 
        %
        % Parameters:
        %  x : N x 1 signal @type Signal
        %  m : M x 1 signal containing the desired (sub)signals @type Signal
        %
        % Return values:
        %  s : M x N selection matrix such that @code m = s*x @endcode
            s = cell2mat(arrayfun(@(xi) ismember(transpose(m),xi),x,'un',0))';
        end
        
        function disp(s)
        % Implements how a signal is displayed in the MATLAB Command Window. 
            fprintf('Signal with UUIDs:\n');
            for k = 1:length(s)
                fprintf('\t%s',num2str([s(k).UUID]));
                if ~isempty(s(k).alias)
                    fprintf(',[%s]\n',num2str([s(k).alias.UUID]));
                else
                    fprintf('\n');
                end
            end
        end
        
        function ids = listaliases(self)
        % Returns a vector with all IDs of the aliases of a signal of
        % length 1.
        %
        % Return values:
        %  ids : vector with all IDs of this signal and its aliases @type double
            assert(length(self)==1,'You can only list all aliases for a signal of length 1.'); 
            a = self.alias;
            ids = self.UUID;
            for i=1:length(a)
                ids = [ids listaliases(a(i))];
            end
            ids = unique(ids);
        end
        
        
        function [sys,conn] = system(self,name)
        % Creates a system realizing a linear combination.
        %
        % Parameters:
        %  name : name of the resulting linear combination model
        %
        % Return values:
        %  sys : the system including the gains and sums required for
        %  obtaining the linear combinations
        %  conn : the connections that are needed corresponding to \c sys
            c = islincomb(self);
            if any(c) && ~all(c)
                error('Cannot make a system of signals of which some are a linear combination and others are not');
            end
            if all(c)
                mod = SSmod(blkdiag(self.multipliers));
                mod.name = name;
                in_ = transpose(horzcat(self.components));
                out_ = self;
                sys = IOSystem(length(in_),length(out_));
                sys = add(sys,mod);
                conn = [in_==sys.in;out_==sys.out];
            else
                sys = [];
                conn = [];
            end
        end
    end
end

