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

classdef NormConstraint < Specification
    % A NormConstraint defines a constraint on a (weighted) norm of a
    % Channel (= a Norm object). The main purpose of this class is to keep
    % track of the upper bound of the norm and to avoid manipulations on
    % the Norm object after the creation of the constraint. 
    
    properties
        norm  = [];   % the norm that is bounded by the constraint @type Norm
        bound = [];   % upper bound on the norm, i.e. @f$\|\cdot\|\leq\text{ bound}@f$ @type double
    end
    
    methods
        function self = NormConstraint(norm,bound)
        % Constructor for NormConstraint objects.
        %
        % Parameters:
        %  norm : the norm that is constrained @type Norm
        %  bound : the upper bound on the norm @type double
        %
        % Return values:
        %  self : the constraint @type NormConstraint
            assert(isnumeric(bound) && isscalar(bound),'The upper bound of a norm should be a numeric scalar.');
            self.norm = norm;
            self.bound = bound;
        end 
        
        function disp(self)
        % Print information
            disp(self.norm);
        end
        
        function le(varargin)
        % Issues an error when the user tries to change the upper bound
        % after creation of the constraint. 
            error('You cannot define an inequality on a constraint.');
        end
        
        function mtimes(varargin)
        % Issues an error when the user tries to reweigh the norm. 
            error('Multiplying a constraint does not make sense.');
        end
        
        function times(varargin)
        % Issues an error when the user tries to reweigh the norm. 
            error('Multiplying a constraint does not make sense.');
        end
        
        function ncarray = plus(n1,n2)
        % Creates a set of constraints that apply to a control problem.
        %
        % Parameters: 
        %  n1 : first constraint @type NormConstraint
        %  n2 : second constraint @type NormConstraint
            ncarray = vertcat(n1,n2);
        end
        
        function ch = ch_p(self)
        % Returns the Channel of the Norm that is constrained. 
        %
        % Return values: 
        %  ch : the channel @type Channel
            ch = self.norm.ch_p; 
        end
        
        function W_out = W_out(self)
        % Returns the input weight of the constrained norm. 
        %
        % Return values: 
        %  W_out : the output weight @type Model
            W_out = self.norm.W_out;
        end
                
        function W_in = W_in(self)
        % Returns the input weight of the constrained norm.
        %
        % Parameters: 
        %  W_in : the input weight @type Model
            W_in = self.norm.W_in; 
        end
        
        function self = setin(self,signal)
        % Sets the intput signals of the channel of the constrained norm.
        %
        % Parameters:
        %   self: the constrained norm @type NormConstraint
        %   signal: the input signals of the norm channel @type Signal
        %
        % Return values:
        %  self : the resulting constrained norm @type NormConstraint
            self.norm = setin(self.norm,signal);
        end
        
        function self = setout(self,signal)
        % Sets the output signals of the channel of the constrained norm.
        %
        % Parameters:
        %   self: the constrained norm @type NormConstraint
        %   signal: the output signals of the norm channel @type Signal
        %
        % Return values:
        %  self : the resulting constrained norm @type NormConstraint
            self.norm = setout(self.norm,signal);
        end
        
        function self = setWout(self,W_out)
        % Sets the intput signals of the channel of the constrained norm.
        %
        % Parameters:
        %   self: the constrained norm @type NormConstraint
        %   W_out: the output weight of the constrained norm @type Model
        %
        % Return values:
        %  self : the resulting constrained norm @type NormConstraint
            self.norm = setWout(self.norm,W_out);
        end
        
        function self = setWin(self,W_in)
        % Sets the output weight of the constrained norm.
        %
        % Parameters:
        %   self: the constrained norm @type NormConstraint
        %   W_in: the input weight of the constrained norm @type Model
        %
        % Return values:
        %  self : the resulting constrained norm @type NormConstraint
            self.norm = setWin(self.norm,W_in);
        end
        
        function self = setchannel(self,ch)
        % Sets the channel of the constrained norm.
        %
        % Parameters:
        %   self: the constrained norm @type NormConstraint
        %   ch: the channel of the constrained norm @type Channel
        %
        % Return values:
        %  self : the resulting constrained norm @type NormConstraint
            self.norm = setchannel(self.norm,ch);
        end
        
        function b = isoutput(self)
        % Checks whether the constrained norm has a non-static output weight.
        %
        % Return values:
        %  b : true if the constrained norm has a non-static output weight @type logical
            b = isoutput(self.norm);
        end
        
        function b = isinput(self)
        % Checks whether the constrained norm has a non-static input weight.
        %
        % Return values:
        %  b : true if the constrained norm has a non-static input weight @type logical
            b = isinput(self.norm);
        end
        
        function b = isstable(self)
        % Checks whether the input and output weights of the constrained
        % Norm are stable or not.
        %
        % Parameters: 
        %  b : true if the weights are all stable @type logical
            b = isstable(self.norm);
        end
        
        function b = isparametric(self)
        % Checks whether the weights of the constrained Norm are parameter-dependent. 
        %
        % Return values:
        %  b : true if the weights are parameter-dependent @type logical
            b = isparametric(self.norm);
        end
        
        function b = isHinf(self)
        % Checks whether the constrained norm is an @f$\mathcal{H}_\infty@f$ norm.
        %
        % Return values:
        %  b : true if the constrained norm is an @f$\mathcal{H}_\infty@f$ norm @type logical
            b = isHinf(self.norm);
        end
        
        function b = isH2(self)
        % Checks whether the constrained norm is an @f$\mathcal{H}_2@f$ norm.
        %
        % Return values:
        %  b : true if the constrained norm is an @f$\mathcal{H}_2@f$ norm @type logical
            b = isH2(self.norm);
        end
        
        function p = normtype(self)
        % Returns the norm type (2 or Inf). 
        %
        % Return values:
        %  p : the norm type (2 or Inf)
            p = normtype(self.norm);
        end
        
        function s = scale(self)
        % Returns the scale of the constrained Norm.  
        %
        % Return values:
        %  s : scale of the constrained Norm @type double
            s = scale(self.norm);
        end
        
        
        function b = upperbound(self)
        % Returns the bound of the constrained Norm.  
        %
        % Return values:
        %  b : upper bound on the constrained Norm @type double
            b = self.bound;
        end
        
        function self = setnorm(self, norm)
        % Sets the norm that is constrained. Only for internal use - end
        % users should not change norms after creation of a constraint.
            self.norm = norm;
        end
        
        function n = getnorm(self)
        % Returns the norm that is constrained.
        %
        % Return values:
        %  n : the norm that is constraint @type Norm
            n = self.norm;
        end
        
        
        function self = setbound(self, bound)
        % Sets the upper bound for the norm. Only for internal use - end
        % users should not change norms after creation of a constraint.
            self.bound = bound;
        end
        
    end

    methods(Static)
        function self = dealscale(self,scale)
        % Includes the scale of the constrained norm in the weights and resets the
        % norm's scale to 1 (or empty). Is a static method since it is also
        % possible to call it on cells. 
        %
        % Parameters:
        %   self: the constraint (\c NormConstraint), or a cell of constraints (\c cell)
        %   scale: optional scale, which the current scale is multiplied
        %   with; 1 if omitted (\c double or \c cell)
        %
        % Return values:
        %  self : the same (set of) constraint(s), but with the scale included in
        %  the weights of the norm (\c NormConstraint or \c cell)
            if iscell(self)
                assert(all(cellfun(@(x) isa(x, 'NormConstraint'), self)), 'dealscale can only operate on NormConstraint objects.');
                normcell = cellfun(@(x) {x.norm}, self);
            else
                normcell = {self.norm};
                self = {self};
            end
            if nargin < 2
                scale = num2cell(ones(size(self)));
            else
                assert(isa(scale,'cell') && length(scale)==length(self), 'scale should be a cell with as many scales as norms.');
            end
            self = cellfun(@(x,y,z) {x.setnorm(Norm.dealscale(y,{z}))}, self, normcell, scale);
        end
    end
end

