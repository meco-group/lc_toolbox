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

classdef Solver
    %SOLVER Abstract solver class
    %   Implements the abstract solver features. Each solver should fill
    %   out the properties of the Solver class as well as the abstract
    %   solve method.
    
    properties
        K;                  % provides the reconstructed controller from the solution.
        H;                  % closed loop model fromed using the controlle robtained.
        performance;        % performance channels.
        
        gamma = [];         % closed loop Hinf upper bound.
        mu = [];            % closed loop H2 upper bound.
        solved = false;     % problem solved or not?.
        
        options;            % options for the solver.
        info;               % captures all residual info coming from the solver.
    end
    
    properties (Access = private,Constant)
        solverlist = {'Hinfsyn','mixedHinfsyn','mixedHinfsynMIMO','mixedHinfsynMIMO_unstab','mixedFixedOrder','systune','SynLPV','mixSynLPV'};
    end
    
    methods   
        function self = setoptions(self,options)
            % Merge defaults and user defined options.
            %
            % Parameters: 
            % options : a structure with user defined options which needs
            % to be merged, for e.g it may contain controller name,
            % explicitly mentioned which solver to use. @type struct
            %
            % Return values:
            % self: options. @type Solver
            
            self.options = mergestruct(options, self.options);
        end
        
        function bodemag(varargin)
            % Overloads MATLAB's Bode magnitude plot, but this one is only 
            % for the desired SISO channels. @type Channel
            %
            % Parameters: 
            % varargin: may contain channels for which bodemag plot is
            % required, and also the additional parameters frequency
            % windows etc. 
            %
            c = cellfun(@(x)isa(x,'Channel'),varargin);
            channels = varargin(c);
            assert(all(cellfun(@issiso,channels)),'Bodemag requires channels to be siso. If you are dealing with mimo responses, use bodemagn');
            solutions = varargin(~c);
            
            if isempty(channels)
                p = cellfun(@(x)horzcat(x.performance.ch_p),solutions,'un',0);
                p = unique(horzcat(p{:}));
                channels = num2cell(p(issiso(p)));
            end
            
            N = length(channels);
            colors = 'bgrcmyk';
            wstyle = '--';

            assert(length(solutions)<=length(colors),'not enough colors to draw');
            colors = colors(1:length(solutions));                
            for k = 1:N
                % Plot the responses
                chan = channels{k};
                m = cellfun(@(x)x.H(chan.out,chan.in),solutions(:),'un',0);
                m = reshape([m';num2cell(colors)],1,[]);

                subplot(1,N,k);
                bodemag(m{:});

                % Plot the weights
                P = bodeoptions('cstprefs'); % check bode preferences for sigma plot
                if strcmp(P.FreqUnits,'Hz'), sc = 2*pi; 
                else sc = 1; end
                lim = axis; % store original axis settings to restore

                [w,b] = cellfun(@(x)x.performance.invsigma(chan.out,chan.in,{sc*1e-10,sc*1e10}),solutions(:),'un',0);
                b = cell2mat(b);
                if any(b)
                    w = w(b);
                    wcolors = cellfun(@(x,y)repmat({[x,wstyle]},1,length(y)),num2cell(colors(b))',w,'un',0);
                    w = reshape([horzcat(w{:});horzcat(wcolors{:})],1,[]); %;horzcat(f{:})
                    hold('on');
                    bodemag(w{:});
                    hold('off');
                    axis(lim); % restore axis settings
                end
                
                title(chan.name);
            end
        end
        
        function sigma(varargin)
            % Overloads MATLAB's Singular value plot (sigma) for desired 
            % channels. @type Channel
            %
            % Parameters: 
            % varargin: may contain channels for which sigma plot is
            % required, and also the additional parameters frequency
            % windows etc. 
            %
            c = cellfun(@(x)isa(x,'Channel'),varargin);
            channels = varargin(c);
            solutions = varargin(~c);
            
            % expand all solutions (might be gridded solutions..)
            solutions = cellfun(@expand,solutions,'un',0);
            solutions = horzcat(solutions{:});
            
            % search for channels if not part of the arguments
            if isempty(channels)
                p = cellfun(@(x)horzcat(x.performance.ch_p),solutions,'un',0);
                p = unique(horzcat(p{:}));
                channels = num2cell(p);
            end
            
            % loop over all channels and solutions
            N = length(channels);
            M = length(solutions);
            colors = 'bgrcmyk';
            
            if length(solutions)>length(colors)
                warning('not enough colors to draw: duplicating colors');
                colors = repmat(colors,[1,ceil(length(solutions)/length(colors))]);
            end
            colors = colors(1:length(solutions));
            
            for k = 1:N
                % Plot the responses
                chan = channels{k};
                m = cell(1,M);
                for j = 1:M
                    s = solutions{j};
                    s.performance = combine(s.performance);
                    Hj = s.H(chan.out,chan.in);
                    c = arrayfun(@(x)eq(x.ch_p,chan),s.performance);
                    if any(c)
                        mj = arrayfun(@(x) cellfun(@(y) x.W_out*y*x.W_in,Hj.content(),'un',0),s.performance(c),'un',0);
                        mj = horzcat(mj{:});
                    else
                        mj = {Hj};
                    end
                    m{j} = reshape([mj;repmat({colors(j)},1,length(mj))],1,[]);
                end
                m = horzcat(m{:});
                
                subplot(1,N,k);
                sigma(m{:});
                title(chan.name);
            end
        end
        
        function list = expand(self)
            % This function is helpful in case of gridded solutions, is
            % used in singular value plot. It creates a list with all
            % possible available solutions.
            % 
            % Return values: 
            % list : cell represting the gridded solution.
            % @type cell
            list = {self};
        end
    end
    
    methods (Static, Abstract)
        capabilities();
    end
    
    methods (Static)
        function [solver] = select(config,specs,vars,options)
            % This function is the most import function in the Solver
            % class, since it is responsible for the automatic and fruitful
            % selection of the solver script for solving the control design problem.
            %
            % Parameters: 
            % config : control configuration @type SystemOfSystems 
            % specs : specifications for controller design. 
            % @type ControllerDesign
            % vars : lists containing system of models. @type cell
            % options : a structure with user defined options which needs
            % to be merged, for e.g it may contain controller name,
            % explicitly mentioned which solver to use. @type struct
            %
            % Return values: 
            % solver : it returns the solver that is selected for the
            % control design problem. @type Solver
            
            assert(length(vars) == 1,'Solving multiple controllers at once not yet implemented');
            % Get rid of the variables and copy the configuration
            [~,r] = vars{1}.empty(); % empty the variables
            config = config.copy(); % copy the config
            restore(r);
            
            % Check for gridded models
            cconf = content(remove_empty(unpack(config)));
            isg = cellfun(@(x)isa(content(x,1),'Gridmod'),cconf);
            if any(isg);
                g = cconf{isg}; gmod = g.content(1);
                g.empty(); g.add(gmod.grid_{1});
                cell_solver = Solver.select(config,specs,vars,options);
                solver = Gridsolver(cell_solver,options);
            else
                solver_list_copy = struct(  'name',{'Hinfsyn','mixedHinfsyn','mixedHinfsynMIMO','mixedHinfsynMIMO_unstab','mixedFixedOrder','systune','HIFOO','SynLPV','mixSynLPV'},...
                                            'constraints',{false,true,true,true,true,true,true,false,true},...
                                            'inout',{0,1,2,2,2,2,2,0,0},...
                                            'norm',{Inf,Inf,Inf,Inf,[2;Inf],[2;Inf],Inf,[2;Inf],[2;Inf]},...
                                            'unstable',{false,false,false,true,false,false,false,false,false},...
                                            'improper',{false,false,false,false,false,false,false,false,false},...
                                            'parametric',{false,false,false,false,false,false,false,true,true},...
                                            'order',{-1,-1,-1,-1,Inf,Inf,Inf,-1,-1});

                % compose solver capabilities
                cap = cellfun(@(x)setfield(eval(['Solver_' x '.capabilities()']),'name',x),Solver.solverlist);

                % Check the solver selection criteria
                criteria.constraints = isCO(specs);
                criteria.inout = isSISO(specs);
                criteria.norm = unique([specs.performance.p]);
                criteria.unstable = ~isstable(specs);
                criteria.improper = false;
                criteria.parametric = isparametric(specs) | any(cellfun(@(x)isparametric(content(x,1)),cconf));

                % Check which solvers can be applied for full order
                remove = (~[cap.constraints]) & criteria.constraints;
                remove = remove | ([cap.inout] < criteria.inout);
                remove = remove | (~arrayfun(@(x) all(ismember(criteria.norm,x.norm)),cap));
                remove = remove | ((~[cap.unstable]) & criteria.unstable);
                remove = remove | ((~[cap.improper]) & criteria.improper);
                remove = remove | ((~[cap.parametric]) & criteria.parametric);
                cap(remove) = [];

                % Choose a solver for the problem
                % Set full order solver
                assert(~isempty(cap),'No appropriate solver found');
                if ~isfield(options,'FullOrderSolver') || ~any(strcmp({cap.name},options.FullOrderSolver))
                    options.FullOrderSolver = cap(1).name; % Choose easiest solution method
                end
                
                % Set the solver object
                if specs.order == -1
                    solver = eval(['Solver_' options.FullOrderSolver '(options)']);
                else
                    % Set reduced order solver
                    cap(~[cap.fixedorder]) = [];
                    assert(~isempty(cap),'No appropriate solver found');
                    if ~isfield(options,'ReducedOrderSolver') || ~any(strcmp({cap.name} ,options.ReducedOrderSolver))
                        options.ReducedOrderSolver = cap(1).name; % Choose easiest solution method
                    end
                    solver = eval(['Solver_' options.ReducedOrderSolver '(options)']);
                end
            end
        end
        
        function [GP,wspecs] = plant(config,specs,vars)
            % This function gather up informations from the control
            % configuration. 
            % 
            % Parameters:
            % config : control configuration. @type SystemOfSystems 
            % specs : specifications for controller design. 
            % @type ControllerDesign
            % vars : lists contining system of models. @type cell
            %
            % Return values: 
            % GP : returns the generalized plant including all exogenous
            % inputs and regulated outputs. @type SystemOfModels
            % wspecs : specifications for controller design to be fed in he
            % solver routine. @type ControllerDesign
            
            % 1. make the openloop -> remove the variables
            [~,r] = cellfun(@empty,vars,'un',0);

            % 2. make systems of the weights
            weights = {};
            connections = {};
            exog_in_ = [];
            exog_out_ = [];
            wspecs = specs;

            % set up weights for augmented plant
            for k = 1:length(specs.performance)
                chan = specs.performance(k).ch_p;
                
                % Add system to generalized plant problem
                if specs.performance(k).isoutput()
                    sysout = IOSystem(specs.performance(k).W_out);
                    weights = [weights,{sysout}];
                    connections = [connections;sysout.in == chan.out];
                    chan.out = sysout.out;
                    wspecs.performance(k).ch_p.out = sysout.out;
                end
                if specs.performance(k).isinput()
                    sysin = IOSystem(specs.performance(k).W_in);
                    weights = [weights,{sysin}];
                    connections = [connections;sysin.out == chan.in];
                    chan.in = sysin.in;
                    wspecs.performance(k).ch_p.in = sysin.in;
                end
                exog_in_ = [exog_in_;chan.in];
                exog_out_ = [exog_out_;chan.out];
            end

            % 3. Throw it all together
            GP = IOSystem(config,weights{:},connections);
            GP = GP.model();
            if ~isparametric(GP.content(1))
                c = ssminreal(simplify(GP.content(1)));
                GP.empty();
                GP.add(c);
            end
            
            exog_in_ = unique(exog_in_);
            exog_out_ = unique(exog_out_);
            GP = GP([exog_out_;specs.ctrl_out],[exog_in_;specs.ctrl_in]);
            
            % 4. Set removed content back
            restore(vertcat(r{:}));
        end
        
        function ch = channels(wspecs)
            % This function is used to set up the specifications of the
            % channels in a structure.
            %
            % Parameters: 
            % wspecs : specifications for controller design to be fed in he
            % solver routine. @type ControllerDesign
            % 
            % Return values:
            % ch : a structure containing specifications of the channels.
            % @type struct
            
            ch.H2 = []; ch.Hinf = []; channel = 0;
            exin = unique(wspecs.in); exout = unique(wspecs.out);
            
            for norm = wspecs.performance                channel = channel + 1;
                if isH2(norm)
                    ch.H2 = [ch.H2 channel];
                elseif isHinf(norm)
                    ch.Hinf = [ch.Hinf channel];
                else
                    error('Unknown type of norm. Not H2 nor Hinf.');
                end
                
                assert(all(ismember(unique(norm.ch_p.in),exin,false)),'Contact the developers');
                assert(all(ismember(unique(norm.ch_p.out),exout,false)),'Contact the developers');
                ch.In{channel} = transpose(selection(unique(norm.ch_p.in),exin));
                ch.Out{channel} = transpose(selection(unique(norm.ch_p.out),exout));
            end
        end
    end
end

