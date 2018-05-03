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

classdef(Abstract) Specification
    % A Specification is an abstract superclass for control problem
    % specifications, that is (for the time being): a Norm (for objective
    % terms), a NormConstraint (as a hard constraint) or an Order (for
    % fixed or reduced-order controller design).
    
    methods
        
        function c = vertcat(varargin)
        % Returns a set of Specification objects as a cell. 
        %
        % Parameters:
        %  varargin : can contain all Specification objects and cells
        %  containing a set of Specification objects 
        %
        % Return values:
        %  c : cell containing all specifications @type cell
        
            isspec = cellfun(@(x) isa(x,'Specification'), varargin);
            iscell = cellfun(@(x) isa(x,'cell'), varargin);
            iscellofspecs = cellfun(@(x) isa(x,'cell') && all(cellfun(@(y) isa(y,'Specification'), x)), varargin);
            rem = ~(isspec | iscell);
            cellclasses = cellfun(@(x) cellfun(@class, x, 'un', 0), varargin((iscell & ~iscellofspecs)), 'un', 0);
            remclass = [cellfun(@class, varargin(rem), 'un', 0) cellclasses{:}];
            
            if ~isempty(remclass) 
                error(['Cannot concatenate specifications or cells of specifications with these types: ' strjoin(remclass)]);
            end
            
            specs = varargin(isspec);
            if ~isempty(varargin(iscellofspecs))
                unpacked_specs = vertcat(varargin{iscellofspecs});
            else
                unpacked_specs = cell(1,0);
            end
            
            c = [specs(:) ; unpacked_specs(:)];
            
        end
        
        function c = horzcat(varargin)
        % Returns a set of Specification objects as a cell.
        %
        % Parameters:
        %  varargin : can contain all Specification objects and cells
        %  containing a set of Specification objects 
        %
        % Return values:
        %  c : cell containing all specifications @type cell
        
            c = vertcat(varargin{:});
        
        end
        
        function c = plus(varargin)
        % Returns a set of Specification objects as a cell.
        %
        % Parameters:
        %  varargin : can contain all Specification objects and cells
        %  containing a set of Specification objects 
        %
        % Return values:
        %  c : cell containing all specifications @type cell
        
            c = vertcat(varargin{:});
        
        end
        
    end
    
end