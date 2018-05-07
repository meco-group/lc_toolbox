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
        
            c = {};
            remclass = {};
            
            for k = 1:length(varargin)
                if isa(varargin{k},'Specification')
                    % stack
                    c = vertcat(c,varargin(k));
                elseif isa(varargin{k},'cell')
                    if all(cellfun(@(x) isa(x,'Specification'),varargin{k}))
                        % unpack and stack
                        unpacked = vertcat(varargin{k}{:});
                        c = [c ; unpacked];
                    else
                        remclass = [remclass ; 'cell'];
                    end
                elseif isempty(varargin{k})
                    % do nothing
                    c = c;
                else
                    % throw an error afterwards
                    remclass = [remclass ; class(varargin{k})];
                end
            end
            
            if ~isempty(remclass)
                error(['Cannot concatenate specifications or cells of specifications with these types: ' strjoin(unique(remclass))]);
            end
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
