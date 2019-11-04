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

classdef (Abstract)ControlProblem < matlab.mixin.Heterogeneous
    %CTRL_PROB
    %   Formulation of the control problem. Norms can be specified both as
    %   objectives and constraints. 
    %   After the computations, sys_K contains the controller, sys_H the
    %   closed loop TF and sys_Z the outputs of the performance channels,
    %   without weighting functions
    
    properties
        % Convenience naming variables
        name              = '';              % Name of the controller
        
        % Performance and plotting channels
        performance       = cell(1,0);       % List of performance criteria to design the controller
        plotting          = cell(1,0);       % List of extra channels which need plotting
        nobj              = 0;               % List of booleans determining the type of performance: objective = true, constraint = false
        
        % Controller K and temporary full order controller
        sys_K             = [];              % Controller transfer function
        ctrl_in           = [];              % List of the plant's control inputs.
        ctrl_out          = [];              % List of the plant's control outputs.
        
        % Closed loop channels with and without weights
        sys_Z                                % Cell array with performance channel's tfs 
        sys_P                                % Closed loop performance of all input/output channels
    
        % Options
        options = struct();
    end
        
    methods
        function obj = ControlProblem(varargin)
            if(nargin>0)
                obj.name = varargin{1};
            end
        end
        
        function in_ = in(self)
            in_ = cellfun(@(x)x.ch_p.in,self.performance,'un',0);
            in_ = vertcat(in_{:});
        end
        
        function out_ = out(self)
            out_ = cellfun(@(x)x.ch_p.out,self.performance,'un',0);
            out_ = vertcat(out_{:});
        end
        
        function nc = nconstr(self)
            nc = length(self.performance) - self.nobj;
        end
        
        %% Control problem formulation
        function obj = addobjective(obj,objective)
            obj = addperformance(obj, objective, true);
        end
        
        function obj = addconstraint(obj,constraint)
            obj = addperformance(obj, constraint, false);
        end
        
        function obj = addplotting(obj, channel)
            obj.plotting(1,size(obj.plotting,2)+(1:length(channel))) = channel;
        end
        
        %% Input/Output classification
        function [isSI,in] = isSingleInput(obj)
            [isSI,in] = isSingleX(obj,'in');
        end
        
        function [isSO,out] = isSingleOutput(obj)
            [isSO,out] = isSingleX(obj,'out');
        end
        
        function [id] = isSISO(obj)
            if((obj.nobj == 1) && (length(obj.performance)==1))
                id = 0;
            else
                [isSI,~] = isSingleInput(obj);
                [isSO,~] = isSingleOutput(obj);
                id = (~isSI) + (~isSO);
            end
        end
        
        function [is_sx,xref] = isSingleX(obj,x)
            for i=1:length(obj.performance)
                if isa(obj.performance{i}, 'Norm')
                    chlist{i} = obj.performance{i}.ch_p.(x);
                elseif isa(obj.performance{i}, 'NormConstraint')
                    chlist{i} = obj.performance{i}.norm.ch_p.(x);
                end
            end
            xref = cellfun(@(x) {unique(x)}, chlist);
            % check equivalence of all inputs
            is_sx = true; j = 2;
            while is_sx && (j <= size(xref,2))
                if length(xref{j}) ~= length(xref{j-1})
                    is_sx = false;
                else
                    is_sx = isalias(xref{j}, xref{j-1});
                end
                j = j+1;
            end
            if is_sx; xref = xref{1}; end
        end
        
        %% Objective/Constraint classification  
        function co = isCO(obj)
            if obj.nobj == length(obj.performance)
                co = false;
            else
                co = true;
            end
        end
        
        function [stable] = isstable(obj)
            stable = all(cellfun(@isstable,obj.performance));
        end
        
        function [isPr] = isProper(obj)
            isPr = true;
            for norm = obj.performance
                if(length(zero(norm.W)) > length(pole(norm.W)))
                    isPr = false;
                    break;
                end
            end
        end
        
        function parametric = isparametric(obj)
            parametric = any(cellfun(@isparametric, obj.performance));
        end
        
        %% Norm type
        function p = norms(obj)
            p1 = cellfun(@(x) x.p, obj.performance(1:obj.nobj));
            p2 = cellfun(@(x) x.norm.p, obj.performance(obj.nobj+1:end));
            p = unique([p1;p2])';
        end
        
        %% Scaling
        function obj = rescale(obj,mode)
            if nargin == 1, mode = 'all'; end
            switch(mode)
                case 'all'
                    obj.performance(1:obj.nobj) = Norm.dealscale(obj.performance(1:obj.nobj));
                case 'obj'
                    objs = 1:obj.nobj;
                    if ~isempty(objs)
                        obj.performance(objs) = Norm.dealscale(obj.performance(objs));
                    end
                case 'constr'
                    constrs = obj.nobj+(1:obj.nconstr());
                    if ~isempty(constrs)
                        obj.performance(constrs) = NormConstraint.dealscale(obj.performance(constrs));
                    end
                case 'constr+bound'
                    constrs = obj.nobj+(1:obj.nconstr());
                    if ~isempty(constrs)
                        obj.performance(constrs) = NormConstraint.dealscale(obj.performance(constrs),cellfun(@(x) {1./x.bound}, obj.performance(constrs)));
                        obj.performance(constrs) = cellfun(@(x) {x.setbound(1)}, obj.performance(constrs));
                    end
            end
        end
    end
    
    methods(Access = private)
        function obj = addperformance(obj,performance,isobjective)
            if isobjective
                obj.performance = [performance ; obj.performance];
                obj.nobj = obj.nobj + length(performance);
            else
                obj.performance = [obj.performance ; performance];
            end
        end
    end
end

