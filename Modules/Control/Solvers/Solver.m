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
    % Solver defines the blueprint for solver interfaces. These are objects
    % that contain all information about the control problem and
    % furthermore have methods available that parse them into an
    % appropriate format for the solver they interface.
    
    properties
        K;                  % the controller @type Model
        H;                  % closed loop system with the obtained controller K @type SystemOfModels
        performance;        % set of all objectives and constraints of type Specification @type cell
        
        gamma = [];         % optimal @f$\mathcal{H}_\infty@f$ norms achieved by the solver @type double
        mu = [];            % optimal @f$\mathcal{H}_2@f$ norms achieved by the solver @type double
        solved = false;     % indicates whether the problem was succesfully solved or not @type logical
        
        options;            % options that are passed to the solver @type struct
        info;               % information that the solver returns @type struct
    end
    
    properties (Access = private,Constant)
        solverlist = {'Hinfsyn','mixedHinfsyn','mixedHinfsynMIMO','mixedHinfsynMIMO_unstab','mixedFixedOrder','systune','HIFOO','SynLPV','mixSynLPV'};
    end
    
    methods   
        function self = setoptions(self,options)
            % Merge default solver settings and user-defined options.
            %
            % Parameters: 
            %  options : a structure with user defined options which needs to be merged @type struct
            %
            % Return values:
            %  self: the solver interface @type Solver
            
            self.options = mergestruct(options, self.options);
        end
        
        function showall(self, varargin)
            % Creates a window with figures containing graphical representations of the
            % objectives, constraints and shaped input-output relations.
            % Only @f$\mathcal{H}_\infty@f$ performance measures are shown
            % (@f$\mathcal{H}_2@f$ norms are hard to represent on a Bode
            % diagram). 
            %
            % Parameters: 
            % varargin: contains Solver objects from different problems 
            
            assert(~verLessThan('matlab','7.8'), 'Visualization of all performance specifications at once is only supported for MATLAB R2009a or higher.');
            colors = 'bgrcmyk';
            allinp = [{self} varargin];
            
            % check inputs
            assert(all(cellfun(@(x) isa(x,'Solver'), allinp)),'The showall call can only handle arguments of type Solver.');
            assert(length(varargin) <= 6, 'Cannot draw more than 7 solutions at once, there are not enough colors for that!');
            
            % check all channels that are in the objects
            allchannels = cellfun(@(x) cellfun(@(y) y.ch_p, x.performance, 'un', 0), allinp, 'un', 0);
            allchannels = vertcat(allchannels{:})'; 
            allchannels = num2cell(unique(Channel.toarray(allchannels)));
            nofch = length(allchannels);
            
            % loop through all channels create a plot of it
            % (strongly inspired by issue 157355 on MATLAB Answers on the MATLAB Central)
             interface = com.mathworks.mde.desk.MLDesktop.getInstance;
             interface.addGroup('allplots');
             interface.setGroupDocked('allplots',0);
             dims = java.awt.Dimension(1,nofch);
             interface.setDocumentArrangement('allplots',2,dims);
             figures = gobjects(1,nofch);
             interface.setDocumentArrangement(['Performance channels: ' strjoin(cellfun(@(x) x.K.name,allinp,'un',0),', ')],1,dims);
             warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
             
            for k=1:nofch
                  thischan = allchannels{k};
                  
                  % preprocess the screen
                   figures(k) = figure('WindowStyle', 'docked', 'Name', sprintf(['Channel: ' thischan.name]), 'NumberTitle', 'off');
                   set(get(handle(figures(k)), 'javaframe'), 'GroupName', ['Performance channels: ' strjoin(cellfun(@(x) x.K.name,allinp,'un',0),', ')]);
                  
                  % the actual plotting
                  l = {}; 
                  for i=1:nargin
                      % new color for this solution
                      thissol = allinp{i};
                      col = colors(i);
                      
                      % check whether the channel occurs in the specs 
                      chs = cellfun(@(x) {getchannel(getnorm(x))}, thissol.performance);
                      chs = cellfun(@(x) x, chs); pl = zeros(size(chs));
                      for j=1:length(chs)
                          if length(chs(j))==length(thischan) && chs(j)==thischan
                              pl(j)=1;
                          else
                              pl(j)=0;
                          end
                      end
                      specs = thissol.performance(logical(pl));
                      
                      if ~isempty(specs)
                          % plot the closed loop performance and the weights
                          if issiso(thischan) % bodemag
                              % the channel itself
                              h = bodeplot(thissol.H(thischan),col);
                              l = [l, {thissol.K.name}];
                              setoptions(h,'PhaseVisible','off');
                              hold on; 

                              % the weights
                              for j=1:length(specs)
                                  if normtype(specs{j})==Inf % 2-norms are not easy to show graphically
                                     bodeplot(inv(specs{j}.W_in*specs{j}.W_out),[col '--']); 
                                     a = gca;
                                     a.Children(1).Annotation.LegendInformation.IconDisplayStyle = 'off'; % don't show constraints in the legend
                                  end
                              end
                          else % sigma
                              l = [l, thissol.K.name];
                              % everything together
                              c = 0;
                              for j=1:length(specs)
                                  if normtype(specs{j})==Inf % 2-norms are not easy to show graphically
                                     sigma(specs{j}.W_out*thissol.H(thischan).content(1)*specs{j}.W_in,col); hold on;
                                     if c>=1
                                        a = gca;
                                        a.Children(1).Annotation.LegendInformation.IconDisplayStyle = 'off'; % avoid same controller multiple times in legend
                                     end
                                     c = c+1;
                                     sigma(ss(1),[col '--']);
                                     a = gca;
                                     a.Children(1).Annotation.LegendInformation.IconDisplayStyle = 'off'; % don't show constraints in the legend
                                  end
                              end
                          end
                      end
                  end
                  
                  % make legend and title
                  drawnow; legend(l,'Location','southeast'); 
                  if issiso(thischan)
                      title(['Bodeplot magnitude: channel ''' thischan.name '''']);
                  else
                      title(['Singular value plot: channel ''' thischan.name '''']);
                  end
            end
            
        end
        
        function bodemag(varargin)
            % Creates a window with figures containing graphical representations of the
            % objectives, constraints and shaped input-output relations.
            % Only @f$\mathcal{H}_\infty@f$ performance measures are shown
            % (@f$\mathcal{H}_2@f$ norms are hard to represent on a Bode
            % diagram). 
            %
            % Parameters: 
            % varargin: contains Solver objects from different problems and
            % Channel objects if only specific channels have to be plotted
            
            assert(~verLessThan('matlab','7.8'), 'Visualization of performance specifications at once is only supported for MATLAB R2009a or higher.');
            colors = 'bgrcmyk';
            
            % check inputs
            assert(all(cellfun(@(x) isa(x,'Solver') || isa(x,'Channel'), varargin)),'The bodemag call can only handle arguments of type Solver and Channel.');
            chs = varargin(cellfun(@(x) isa(x,'Channel'), varargin)); 
            sols = varargin(cellfun(@(x) isa(x,'Solver'), varargin));
            assert(length(sols) <= 6, 'Cannot draw more than 7 solutions at once, there are not enough colors for that!');
            
            % check all channels that are in the objects
            allchs = cellfun(@(x) cellfun(@(y) y.ch_p, x.performance, 'un', 0), sols, 'un', 0);
            allchs = vertcat(allchs{:})'; 
            allchs = unique(Channel.toarray(allchs));
            chs = cellfun(@(x) x, chs);
            if isempty(chs)
                chs = num2cell(allchs);
            else
                [lia,locb] = ismember(chs,allchs);
                assert(any(lia), 'You requested to draw a channel that does not belong to any of the solutions you provided.');
                chs = num2cell(allchs(locb));
            end
            idxsiso = cellfun(@issiso,chs);
            if ~all(idxsiso)
                warning('Only SISO channels are drawn on a bodemag performance plot. MIMO channels are ignored; use sigma for them instead.');
            end
            chs = chs(idxsiso);
            nofch = length(chs);
            
            % loop through all channels create a plot of it
            % (strongly inspired by issue 157355 on MATLAB Answers on the MATLAB Central)
             interface = com.mathworks.mde.desk.MLDesktop.getInstance;
             interface.addGroup('allplots');
             interface.setGroupDocked('allplots',0);
             dims = java.awt.Dimension(1,nofch);
             interface.setDocumentArrangement('allplots',2,dims);
             figures = gobjects(1,nofch);
             interface.setDocumentArrangement(['Performance channels: ' strjoin(cellfun(@(x) x.K.name,sols,'un',0),', ')],1,dims);
             warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
             
            for k=1:nofch
                  thischan = chs{k};
                  
                  % preprocess the screen
                   figures(k) = figure('WindowStyle', 'docked', 'Name', sprintf(['Channel: ' thischan.name]), 'NumberTitle', 'off');
                   set(get(handle(figures(k)), 'javaframe'), 'GroupName', ['Performance channels: ' strjoin(cellfun(@(x) x.K.name,sols,'un',0),', ')]);
                  
                  % the actual plotting
                  l = {}; 
                  for i=1:length(sols)
                      % new color for this solution
                      thissol = sols{i};
                      col = colors(i);
                      
                      % check whether the channel occurs in the specs 
                      chs_ = cellfun(@(x) {getchannel(getnorm(x))}, thissol.performance);
                      chs_ = cellfun(@(x) x, chs_); pl = zeros(size(chs_));
                      for j=1:length(chs_)
                          if length(chs_(j))==length(thischan) && chs_(j)==thischan
                              pl(j)=1;
                          else
                              pl(j)=0;
                          end
                      end
                      specs = thissol.performance(logical(pl));
                      
                      if ~isempty(specs)
                          % plot the closed loop performance and the weights
                              % the channel itself
                              h = bodeplot(thissol.H(thischan),col);
                              l = [l, {thissol.K.name}];
                              setoptions(h,'PhaseVisible','off');
                              hold on; 

                              % the weights
                              for j=1:length(specs)
                                  if normtype(specs{j})==Inf % 2-norms are not easy to show graphically
                                     bodeplot(inv(specs{j}.W_in*specs{j}.W_out),[col '--']); 
                                     a = gca;
                                     a.Children(1).Annotation.LegendInformation.IconDisplayStyle = 'off'; % don't show constraints in the legend
                                  end
                              end
                      end
                  end
                  
                  % make legend and title
                  drawnow; legend(l,'Location','southeast'); 
                  title(['Bodeplot magnitude: channel ''' thischan.name '''']);
            end
        end
        
        function sigma(varargin)
            % Creates a window with figures containing graphical representations of the
            % objectives, constraints and shaped input-output relations.
            % Only @f$\mathcal{H}_\infty@f$ performance measures are shown
            % (@f$\mathcal{H}_2@f$ norms are hard to represent on a
            % singular value plot). 
            %
            % Parameters: 
            % varargin: contains Solver objects from different problems and
            % Channel objects if only specific channels have to be plotted
            
            assert(~verLessThan('matlab','7.8'), 'Visualization of performance specifications at once is only supported for MATLAB R2009a or higher.');
            colors = 'bgrcmyk';
            
            % check inputs
            assert(all(cellfun(@(x) isa(x,'Solver') || isa(x,'Channel'), varargin)),'The bodemag call can only handle arguments of type Solver and Channel.');
            chs = varargin(cellfun(@(x) isa(x,'Channel'), varargin)); 
            sols = varargin(cellfun(@(x) isa(x,'Solver'), varargin));
            assert(length(sols) <= 6, 'Cannot draw more than 7 solutions at once, there are not enough colors for that!');
            
            % check all channels that are in the objects
            allchs = cellfun(@(x) cellfun(@(y) y.ch_p, x.performance, 'un', 0), sols, 'un', 0);
            allchs = vertcat(allchs{:})'; 
            allchs = unique(Channel.toarray(allchs));
            chs = [chs{:}]';
            if isempty(chs)
                chs = num2cell(allchs);
            else
                [lia,locb] = ismember(chs,allchs);
                assert(all(lia), 'You requested to draw a channel that does not belong to any of the solutions you provided.');
                chs = num2cell(allchs(locb));
            end
            nofch = length(chs);
            
            % loop through all channels create a plot of it
            % (strongly inspired by issue 157355 on MATLAB Answers on the MATLAB Central)
             interface = com.mathworks.mde.desk.MLDesktop.getInstance;
             interface.addGroup('allplots');
             interface.setGroupDocked('allplots',0);
             dims = java.awt.Dimension(1,nofch);
             interface.setDocumentArrangement('allplots',2,dims);
             figures = gobjects(1,nofch);
             interface.setDocumentArrangement(['Performance channels: ' strjoin(cellfun(@(x) x.K.name,sols,'un',0),', ')],1,dims);
             warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
             
            for k=1:nofch
                  thischan = chs{k};
                  
                  % preprocess the screen
                   figures(k) = figure('WindowStyle', 'docked', 'Name', sprintf(['Channel: ' thischan.name]), 'NumberTitle', 'off');
                   set(get(handle(figures(k)), 'javaframe'), 'GroupName', ['Performance channels: ' strjoin(cellfun(@(x) x.K.name,sols,'un',0),', ')]);
                  
                  % the actual plotting
                  l = {}; 
                  for i=1:length(sols)
                      % new color for this solution
                      thissol = sols{i};
                      col = colors(i);
                      
                      % check whether the channel occurs in the specs 
                      chs_ = cellfun(@(x) {getchannel(getnorm(x))}, thissol.performance);
                      chs_ = cellfun(@(x) x, chs); pl = zeros(size(chs_));
                      for j=1:length(chs_)
                          if length(chs_(j))==length(thischan) && chs_(j)==thischan
                              pl(j)=1;
                          else
                              pl(j)=0;
                          end
                      end
                      specs = thissol.performance(logical(pl));
                      
                      if ~isempty(specs)
                          % plot the closed loop performance and the weights
                          l = [l, thissol.K.name];
                          % everything together
                          c = 0;
                          for j=1:length(specs)
                              if normtype(specs{j})==Inf % 2-norms are not easy to show graphically
                                 sigma(specs{j}.W_out*thissol.H(thischan).content(1)*specs{j}.W_in,col); hold on;
                                 if c>=1
                                     a = gca;
                                     a.Children(1).Annotation.LegendInformation.IconDisplayStyle = 'off'; % avoid same controller multiple times in legend
                                 end
                                 c = c+1;
                                 sigma(ss(1),[col '--']);
                                 a = gca;
                                 a.Children(1).Annotation.LegendInformation.IconDisplayStyle = 'off'; % don't show constraints in the legend
                              end
                          end
                      end
                  end
                  
                  % make legend and title
                  drawnow; legend(l,'Location','southeast'); 
                  title(['Singular value plot: channel ''' thischan.name '''']);
            end
        end
    end
    
    methods (Static, Abstract)
        capabilities();
    end
    
    methods (Static)
        function [solver] = select(config,specs,vars,options)
            % Selects an appropriate solver amongst the interfaced solvers
            % to solve the control problem. 
            %
            % Parameters: 
            %  config : control configuration @type SystemOfSystems 
            %  specs : specifications for controller design @type ControllerDesign
            %  vars : lists containing all SystemOfModels objects that are variables in the problem (i.e. the controller(s)) @type cell
            %  options : options the user wants to pass to the solver and
            %  to the solver selection routine (e.g. to enforce the use of a certain
            %  solver) @type struct
            %
            % Return values: 
            %  solver : returns a Solver interface for the selected solver @type Solver
            
            assert(length(vars) == 1,'Solving multiple controllers at once not yet implemented');
            
            % Get rid of the variables and copy the configuration
            [~,r] = vars{1}.empty(); % empty the variables
            config = config.copy(); % copy the config
            restore(r);
            
            % Check for gridded models
            cconf = content(remove_empty(unpack(config)));
            isg = cellfun(@(x)isa(content(x,1),'Gridmod'),cconf);
            if any(isg)
                g = cconf{isg}; gmod = g.content(1);
                g.empty(); g.add(gmod.grid_{1});
                cell_solver = Solver.select(config,specs,vars,options);
                solver = Gridsolver(cell_solver,options);
            else
                
                % Compose solver capabilities (loop through all solvers of
                % 'solverlist' and request its capabilities)
                cap = cellfun(@(x)setfield(eval(['Solver_' x '.capabilities()']),'name',x),Solver.solverlist);

                % Check the solver selection criteria
                [~,Woutus,Winus,ch] = Solver.plant(config,specs,vars);
                criteria.constraints = isCO(specs);
                criteria.inout = ~isequal(ch.In{1},ch.In{:}) + ~isequal(ch.Out{1},ch.Out{:});
                criteria.norm = norms(specs);
                criteria.unstable = ~(isempty(Woutus.a) && isempty(Winus.a));
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
                assert(~isempty(cap),'No appropriate solver found');
                if ~isfield(options,'FullOrderSolver') || ~any(strcmp({cap.name},options.FullOrderSolver))
                    options.FullOrderSolver = cap(1).name; % Choose easiest solution method
                end
                
                % Set the solver object
                if specs.order == -1
                    % Set full order solver
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
        
        function [GP,Woutus,Winus,ch] = plant(config,specs,vars)
            % Constructs the generalized plant based on the control 
            % configuration the user specified and on the specifications.
            % 
            % Parameters:
            %  config : the control configuration @type SystemOfSystems 
            %  specs : specifications of the control problem @type ControllerDesign
            %  vars : list containing the variables (controllers that are sought for, SystemOfModels) @type cell
            %  minimal : if true, multiple occurences of identical exogenous inputs or performance outputs are avoided in the generalized plant @type logical
            %
            % Return values: 
            %  GP : returns the generalized plant including all exogenous inputs and regulated outputs. @type SystemOfModels
            %  wspecs : adapted specifications, such that the input signal of every channel becomes the input signal of the corresponding input weight and the output signal of every signal becomes the output of its corresponding output weight @type ControllerDesign
                        
            % 1. Make the openloop -> remove the variables
            [~,r] = cellfun(@empty,vars,'un',0);

            % 2. Make systems of the weights
            norms = specs.performance;
            cellfun(@disp,norms);

            % check which weighted inputs are equal
            [Win,connections_in,norms] = Norm.tofilter(norms,'in');
            [Wout,connections_out,norms] = Norm.tofilter(norms,'out');
            
            % make channel structures
            cellfun(@disp,norms);
            ch = Solver.channels(norms,Wout,Win);
            
            % 3. Throw it all together
            GP = config([vertcat(connections_out{:,2});specs.ctrl_out],[vertcat(connections_in{:,2});specs.ctrl_in]);
            GP = GP.content(1);
            
            % put filter splitting and multiplication here
            try
                [Wins,Winus] = stabsep(std(Win.content(1)));
                if ~isempty(Winus.a) % unstable modes in input filter
                    Winus = [Winus;ss(eye(size(Win,2)))];
                    Wins = [ss(eye(size(Win,1))),Wins];
                else
                    Wins = std(Win.content(1));
                    Winus = ss(eye(size(Wins,2)));
                end
            catch
                warning('stabsep cannot separate modes for improper systems. The stable/unstable decomposition is omitted for now.');
                Wins = std(Win.content(1));
                Winus = ss(eye(size(Wins,2)));
            end
              
            try
                [Wouts,Woutus] = stabsep(std(Wout.content(1)));
                if ~isempty(Woutus.a) % unstable modes in output filter
                    Wouts = [ss(eye(size(Wout,2)));Wouts];
                    Woutus = [Woutus,ss(eye(size(Wout,1)))];
                else
                    Wouts = std(Wout.content(1));
                    Woutus = ss(eye(size(Wouts,1)));
                end
            catch
                warning('stabsep cannot separate modes for improper systems. The stable/unstable decomposition is omitted for now.');
                Wouts = std(Wout.content(1));
                Woutus = ss(eye(size(Wouts,1)));
            end
            
            % multiply the stable parts
            Wouts_ = blkdiag(Wouts,ss(eye(length(specs.ctrl_out))));
            Wins_ = blkdiag(Wins,ss(eye(length(specs.ctrl_in))));
            GP = fromstd(Wouts_)*GP*fromstd(Wins_);
            if ~isparametric(GP)
                GP = ssminreal(simplify(GP));
            end
            
            % 4. Set removed content back
            restore(vertcat(r{:}));
        end
        
        function ch = channels(norms,Wout,Win)
            % Returns a structure defining the performance channels.
            %
            % Parameters: 
            %  wspecs : specifications for the controller design, where the norm inputs and outputs are redefined as in the generalized plant(= output of Solver::plant) @type ControllerDesign
            % 
            % Return values:
            %  ch : structure defining the performance channels channels @type struct
            
            ch.H2 = []; ch.Hinf = []; channel = 0;
            exin = Win.in; exout = Wout.out;
            for k = 1:length(norms)
                spec = norms{k};
                channel = channel + 1;
                if isH2(spec)
                    ch.H2 = [ch.H2 channel];
                elseif isHinf(spec)
                    ch.Hinf = [ch.Hinf channel];
                else
                    error('Unknown type of norm. Not H2 nor Hinf.');
                end
                
                ch_p = spec.ch_p;
                assert(all(ismember(unique(ch_p.in),exin,false)),'Contact the developers');
                assert(all(ismember(unique(ch_p.out),exout,false)),'Contact the developers');
                ch.Out{channel} = transpose(selection(unique(ch_p.out),exout));
                ch.In{channel} = selection(unique(ch_p.in),exin);
            end
        end
    end
end

