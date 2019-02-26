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

classdef (InferiorClasses = {?zpk,?tf,?ss,?frd}) AbstractDSSmod < AbstractLFTmod
% AbstractDSSmod provides a blueprint consisting of the properties and
% methods every DSS model of LCToolbox should have. 
        
    methods
        function self = AbstractDSSmod(A,B,C,D,E,varargin)
        % Constructor for AbstractDSSmod objects. 
        %
        % Parameters: 
        %  A : the system's A matrix, either \c double or \c Function (in
        %  case of parameter dependency)
        %  B : the system's B matrix, either \c double or \c Function (in
        %  case of parameter dependency)
        %  C : the system's C matrix, either \c double or \c Function (in
        %  case of parameter dependency)
        %  D : the system's D matrix, either \c double or \c Function (in
        %  case of parameter dependency)
        %  E : the system's D matrix, either \c double or \c Function (in
        %  case of parameter dependency) 
        %  varargin : can contain additional information that is passed to AbstractLFTmod::AbstractLFTmod
        
            assert(all(size(A,1) == [size(A,2),size(B,1),size(C,2)]),'state dimensions inconsistent');
            assert(size(B,2)==size(D,2),'input dimensions inconsistent');
            assert(size(C,1)==size(D,1),'output dimensions inconsistent');
            if isempty(E), E = eye(size(A)); end
            
            % construct M
            nx = size(A,1); nu = size(B,2); ny = size(C,1);
            M = {[],zeros(0,nx),zeros(0,nu),[];...
                 zeros(nx,0),A,B,zeros(nx,0);...
                 zeros(ny,0),C,D,zeros(ny,0);...
                 [],zeros(0,nx),zeros(0,nu),[]};
            self@AbstractLFTmod(M,[],[],E,varargin{:});
            check(self);
        end
        
        function [A] = A(self,varargin)
        % Returns the model's A matrix.
        %
        % Parameters:
        %  varargin : allows subscripting of the A matrix 
        %
        % Return values: 
        %  A : the model's A matrix, either \c double or \c Function
            A = self.M{2,2};
            if nargin>1, A = A(varargin{:}); end
        end
        
        function [B] = B(self,varargin)
        % Returns the model's B matrix.
        %
        % Parameters:
        %  varargin : allows subscripting of the B matrix 
        %
        % Return values: 
        %  B : the model's A matrix, either \c double or \c Function
            B = self.M{2,3};
            if nargin>1, B = B(varargin{:}); end
        end
        
        function [C] = C(self,varargin)
        % Returns the model's C matrix.
        %
        % Parameters:
        %  varargin : allows subscripting of the C matrix 
        %
        % Return values: 
        %  C : the model's C matrix, either \c double or \c Function
            C = self.M{3,2};
            if nargin>1, C = C(varargin{:}); end
        end
        
        function [D] = D(self,varargin)
        % Returns the model's D matrix.
        %
        % Parameters:
        %  varargin : allows subscripting of the D matrix 
        %
        % Return values: 
        %  D : the model's D matrix, either \c double or \c Function
            D = self.M{3,3};
            if nargin>1, D = D(varargin{:}); end
        end
        
        function [a] = a(self,varargin)
        % Returns the model's A matrix. 
        % Lower case version of AbstractDSSmod::A.
            a = A(self,varargin{:});
        end
        
        function [b] = b(self,varargin)
        % Returns the model's B matrix. 
        % Lower case version of AbstractDSSmod::B.
            b = B(self,varargin{:});
        end
        
        function [c] = c(self,varargin)
        % Returns the model's C matrix. 
        % Lower case version of AbstractDSSmod::C.
            c = C(self,varargin{:});
        end
        
        function [d] = d(self,varargin)
        % Returns the model's D matrix. 
        % Lower case version of AbstractDSSmod::D.
            d = D(self,varargin{:});
        end
        
        function [A,B,C,D,E] = getdssdata(self)
        % Returns the model's descriptor state-space matrices.
        %
        % Return values:
        %  A : the model's A matrix, either \c double or \c Function
        %  B : the model's B matrix, either \c double or \c Function
        %  C : the model's C matrix, either \c double or \c Function
        %  D : the model's D matrix, either \c double or \c Function
        %  E : the model's E matrix, either \c double or \c Function
            A = self.A;
            B = self.B;
            C = self.C;
            D = self.D;
            E = self.E;
        end
        
        function self = setdssdata(self,A,B,C,D,E)
        % Sets the model's descriptor state-space matrices.
        %
        % Parameters:
        %  A : the model's A matrix, either \c double or \c Function (in
        %  case of parameter dependency)
        %  B : the model's B matrix, either \c double or \c Function (in
        %  case of parameter dependency)
        %  C : the model's C matrix, either \c double or \c Function (in
        %  case of parameter dependency)
        %  D : the model's D matrix, either \c double or \c Function (in
        %  case of parameter dependency)
        %  E : the model's E matrix, either \c double or \c Function (in
        %  case of parameter dependency)
        %
        % Return values:
        %  self : returns the same object with new descriptor state-space matrices @type AbstractDSSmod
            nx = size(A,1); nu = size(B,2); ny = size(C,1);
            M = {[],zeros(0,nx),zeros(0,nu),[];...
                 zeros(nx,0),A,B,zeros(nx,0);...
                 zeros(ny,0),C,D,zeros(ny,0);...
                 [],zeros(0,nx),zeros(0,nu),[]};
            
            self = setdata(self,M,[],[],E);
            check(self);
        end
        
        function [As,Bs,Cs,Ds,Es] = grid_eval(self,grid,args)
        % Evaluates the model for a grid of parameter values.
        %
        % Parameters:
        %  grid : cell with grid{i} a vector (\c double) with the gridpoints for
        %  scheduling parameter i
        %  args : cell with args{i} the name (\c char) of scheduling parameter i @type cell
        %
        % Return values:
        %  As : array of A matrices with the appropriate size @type double
        %  Bs : array of B matrices with the appropriate size @type double
        %  Cs : array of C matrices with the appropriate size @type double
        %  Ds : array of D matrices with the appropriate size @type double
        %  Es : array of E matrices with the appropriate size @type double
            nargs = length(args);
            As = permute(safegrideval(self.A,grid,args),[nargs+(1:2),1:nargs]);
            Bs = permute(safegrideval(self.B,grid,args),[nargs+(1:2),1:nargs]);
            Cs = permute(safegrideval(self.C,grid,args),[nargs+(1:2),1:nargs]);
            Ds = permute(safegrideval(self.D,grid,args),[nargs+(1:2),1:nargs]);
            Es = permute(safegrideval(self.E,grid,args),[nargs+(1:2),1:nargs]);
        end
        
        function [As,Bs,Cs,Ds,Es] = eval(self,val,args)
        % Evaluates the model for one combination of scheduling parameters.
        %
        % Parameters:
        %  val : vector with val(i) the value of scheduling parameter i
        %  @type double
        %  args : cell with args{i} the name (\c char) of scheduling parameter i @type cell
        %
        % Return values:
        %  As : A matrix for the requested parameter combination @type double
        %  Bs : B matrix for the requested parameter combination @type double
        %  Cs : C matrix for the requested parameter combination @type double
        %  Ds : D matrix for the requested parameter combination @type double
        %  Es : E matrix for the requested parameter combination @type double
            if size(val,1) == length(val); val = val'; end
            if ~iscell(val); grid = num2cell(val); else; grid = val; end
            [As,Bs,Cs,Ds,Es] = grid_eval(self,grid,args);
        end
        
    end
end

