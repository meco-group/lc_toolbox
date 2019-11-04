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

classdef AbstractLPVmod 
% AbstractLPVmod provides a blueprint consisting of the properties and
% methods every LPV model of LCToolbox should have. 
    
    properties (Access=private)
        parameters_ = {};   % cell containing the scheduling parameters of the model @type cell
    end
        
    methods
        function self = AbstractLPVmod(parameters)
        % Constructor for AbstractLPVmod objects. 
        %
        % Parameters: 
        %  parameters : cell with scheduling parameters (\c SchedulingParameter objects) that are in the model @type cell
            if ~iscell(parameters)
                parameters = {parameters};
            end
            self.parameters_ = parameters;
        end
        
        function b = isstable(varargin)
        % Checks the stability of the LPV model.
        % Overloads MATLAB's \c isstable() for AbstractLPVmod objects.
            varargin = stdargs(varargin);
            b = isstable(varargin{:});
        end
        
        function b = isparametric(self)
        % Checks whether the model is parameter dependent. 
        %
        % Return values: 
        %  b : boolean reflecting whether the model is parameter dependent
        %  @type logical
            b = (self.nparameters > 0);
        end
        
        function p = parameters(self)
        % Returns a cell with the scheduling parameters of the model. 
        %
        % Return values: 
        %  p : cell with the scheduling parameters of the model @type cell
            p = self.parameters_;
        end
        
        function a = arguments(self)
        % Returns a cell with the names of the scheduling parameters of the model. 
        %
        % @note This is method that is required in order to ease the interface 
        % with the Optispline toolbox. 
        %
        % Return values: 
        %  p : cell with the names (\c char) of the scheduling parameters @type cell
            a = cellfun(@(x) x.tensor_basis.arguments, self.parameters,'un',0);
        end
        
        function n = nparameters(self)
        % Returns the number of scheduling parameters of the model.
        %
        % Return values: 
        %  n : number of scheduling parameters of the model @type double
            n = length(self.parameters_);
        end
        
        function p = mergeparameters(varargin)
        % Removes duplicates in the list of scheduling parameters of the
        % model.
        %
        % Return values: 
        %  p : list of all scheduling parameters of the model without
        %  duplicates @type cell
            p = cellfun(@parameters,varargin,'un',0);
            p = horzcat(p{:});
            [~,ip,~] = unique(cellfun(@(x)arguments(tensor_basis(x)),p,'un',0));
            p = p(ip);
        end
                
        function [grid,args] = makegrid(self,varargin) 
        % Grids the scheduling parameter space.
        %
        % Parameters :
        % varargin : supported calls for this function are:
        % - for the univariate case (only 1 parameter):
        %   - @verbatim makegrid(self) @endverbatim creates a default equidistant grid (n = 5)
        %   - @verbatim makegrid(self,n) @endverbatim creates an equidistant grid with n points
        %   - @verbatim makegrid(self,[g]) @endverbatim uses the vector g (\c double) of gridpoints that is provided
        %   - @verbatim makegrid(self,g) @endverbatim same as previous call with [g] = g{1,1}
        % - for the multivariate case (multiple parameters):
        %   - @verbatim makegrid(self) @endverbatim creates a default equidistant grid (n = 5) in each dimension
        %   - @verbatim makegrid(self,n) @endverbatim creates an equidistant grid with n points along each dimension
        %   - @verbatim makegrid(self,v,args) @endverbatim uses v (\c cell)
        %   representing a grid such that v{i} is a vector (\c double) with
        %   the gridpoints for scheduling parameters with names args{i}
        %
        % Return values: 
        %  grid : grid with the values of the parameters @type double
        %  args : names of the scheduling parameters, in the same order as \c grid @type char
            
            if nargin == 1, varargin{1} = 5; end %std density
            
            % arguments %
            if nargin < 3
                args = arguments(self); 
            else
                args = varargin{2};
            end
            
            % grid %
            if iscell(varargin{1})
                grid = varargin{1};
            else
                if nparameters(self) < length(varargin{1})
                    grid = varargin(1);
                else
                    if nparameters(self) ~= length(varargin{1})
                        varargin{1} = repmat(varargin{1},[1,length(self.parameters())]);
                    end
                    grid = parameter_grid(self,varargin{1});
                end
            end
        end
        
         % tools
%         function varargout = lsim(varargin) % for the time being not used
%         as default implementation
%             [varargout{1:nargout}] = lsimlpv(varargin{:});
%         end
        
        function [grid] = parameter_grid(self,density)
        % Creates a linearly spaced grid by sampling the parameter space.
        % Each dimension can have a different density. 
        %
        % Parameters:
        %  density : vector containing the number of gridpoints for every parameter (same order as cell of \c parameters_) @type double
        %
        % Return values: 
        %  grid : cell of which each column contains the sample values of its corresponding parameter @type cell
            assert(length(density) == nparameters(self),'parameter_grid expects the grid density to have as many entries as parameter');
            grid = cell(1,self.nparameters());
            for k = 1:self.nparameters()
                grid{1,k} = linspace(self.parameters{k}.basis.domain.min+sqrt(eps),self.parameters{k}.basis.domain.max-sqrt(eps),density(k));
            end
        end
        
    end
end

