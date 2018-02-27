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
        performance       = Norm.empty(1,0); % List of performance criteria to design the controller
        plotting          = Norm.empty(1,0); % List of extra channels which need plotting
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
            in_ = arrayfun(@(x)x.ch_p.in,self.performance,'un',0);
            in_ = vertcat(in_{:});
        end
        
        function out_ = out(self)
            out_ = arrayfun(@(x)x.ch_p.out,self.performance,'un',0);
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
            xref = unique(obj.performance(1,1).ch_p.(x));
            % check equivalence of all inputs
            is_sx = true; j = 2;
            while is_sx && (j <= size(obj.performance,2))
                temp = unique(obj.performance(1,j).ch_p.(x));
                if length(xref) ~= length(temp)
                    is_sx = false;
                else
                    is_sx = isalias(temp, xref);
                end
                j = j+1;
            end
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
            stable = all(arrayfun(@isstable,obj.performance));
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
            parametric = false;
            for norm = obj.performance
                parametric = parametric | isparametric(norm);
            end
        end
        
        %% Compute generalized plant
        function [gp_prob_aug] = constructGeneralizedPlant(obj, gp_prob)
%             % Not a clean copy
%             gp_prob_aug = gp_prob;
%             gp_prob_aug.exog_out = []; gp_prob_aug.exog_in = [];
%             gp_prob_aug.ctrl_out = obj.ctrl_out; gp_prob_aug.ctrl_in = obj.ctrl_in; 
%             
%             % set up weights for augmented plant
%             for k = 1:length(obj.performance)
%                 % Include scale in solution
%                 obj.performance(k) = unify(obj.performance(k));
%                 
%                 % Add system to generalized plant problem
%                 if obj.performance(k).isoutput()
%                     in = obj.performance(k).TF.out;
%                     out = Signal(size(obj.performance(k).W,1));
%                     gp_prob_aug = gp_prob_aug.setz([gp_prob_aug.exog_out;out]);
%                     gp_prob_aug = gp_prob_aug.setw(union(gp_prob_aug.exog_in,obj.performance(k).TF.in));
%                 else
%                     in = Signal(size(obj.performance(k).W,2));
%                     out = obj.performance(k).TF.in;
%                     gp_prob_aug = gp_prob_aug.setw([gp_prob_aug.exog_in;in]);
%                     gp_prob_aug = gp_prob_aug.setz(union(gp_prob_aug.exog_out,obj.performance(k).TF.out));
%                 end
%                 gp_prob_aug = addsystem(gp_prob_aug,fromstd(obj.performance(k).W,in,out));
%             end
            

        end
        
        
        
        %% Close the loop
        function obj = closeLoop(obj,gp_prob_uw,norms)
            % Insert the controller as system with input the control output
            % and as output the control input of the plant
            obj.sys_K = setio(obj.sys_K,obj.ctrl_out,obj.ctrl_in);
            gp_prob_uw = addsystem(gp_prob_uw,obj.sys_K);
            
            % Remove all assigned inputs/outputs for automatic assignment
            gp_prob_uw.exog_in = [];
            gp_prob_uw.exog_out = [];
            gp_prob_uw.ctrl_in = [];
            gp_prob_uw.ctrl_out = [];

            % Build the plant and make an LTIsys model out of it
            gp_prob_uw = build(gp_prob_uw);
%             obj.sys_P = LTIsys(gp_prob_uw.sys_G,gp_prob_uw.exog_in,gp_prob_uw.exog_out);
            obj.sys_P = gp_prob_uw.sys_G;

            % Retrieve the performance tfs
            % TODO: fix implicit linkage of norms and resulting tfs
            obj.sys_Z = cell(length(norms),1);
            for k = 1:length(norms)
                obj.sys_Z{k,1} = obj.sys_P(unique(norms(1,k).TF.out,'stable'),unique(norms(1,k).TF.in,'stable'));
            end
        end
        
        %% Channels
        function [channels] = getChannels(obj,varargin)
            % Set up norms and channel variables
            list = [obj.performance obj.plotting];
            channels = struct('io',{cell(1,2)},'name',{},'frf',{},'weights',{{}});
            
            id = 0; %Channel id -> contains number of channels in the end
            k = 0;
            
            for norms = list 
                k = k+1; i = 0;
                [normsplit,frfs] = unstack(norms,obj.sys_Z{k,1});
                
                for norm = normsplit
                    match = 0; l = 1;                
                    while((l<=id)&&(~match)) 
                        % Check for transfer function equality - take
                        % equivalence into account
                        if(nargin>1)
                            [match,~] = abstractEquality(norm.TF,Channel(channels{l}.io{1},channels{l}.io{2}),varargin{1});
                        else
                            match = (norm.TF == Channel(channels{l}.io{1},channels{l}.io{2}));
                        end
                        l = l+1;
                    end
                    if(match)
                        % Add weight to channel if there is a match: no new
                        % channel found
                        channels{l-1}.weights = [channels{l-1}.weights;{norm.W}];
                    else
                        % Add a new channel to the list: new channel discovered
                        id = id+1; i = i+1;
                        channels{id}.io = {norm.TF.in norm.TF.out};
                        channels{id}.name = norm.TF.name;
                        channels{id}.frf = frfs{i,1};
                        channels{id}.weights = {norm.W};
                    end
                end
            end
            
            % remove all empty channels
            for k = 1:id
                for j = length(channels{k}.weights):-1:1
                    if(isempty(channels{k}.weights{j}))
                        channels{k}.weights(j) = [];
                    end
                end                
            end
        end
        
        function obj = rescale(obj,mode)
            if nargin == 1, mode = 'all'; end
            switch(mode)
                case 'all'
                    obj.performance = dealscale(obj.performance);
                case 'obj'
                    objs = 1:obj.nobj;
                    if ~isempty(objs)
                        obj.performance(objs) = dealscale(obj.performance(objs));
                    end
                case 'constr'
                    constrs = obj.nobj+(1:obj.nconstr());
                    if ~isempty(constrs)
                        obj.performance(constrs) = dealscale(obj.performance(constrs));
                    end
            end
        end
    end
    
    methods(Access = private)
        function obj = addperformance(obj,performance,isobjective)
            if isobjective
                obj.performance = [performance obj.performance];
                obj.nobj = obj.nobj + length(performance);
            else
                obj.performance = [obj.performance performance];
            end
        end
    end
end

