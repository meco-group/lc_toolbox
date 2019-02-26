classdef Solver_mixedpoleMIMO < Solver
    methods
        function self = Solver_mixedpoleMIMO(options)
            disp('Solving with mixedpoleMIMO, because there is a closed loop pole regions constraint added.')
            
            if nargin > 0
                self = setoptions(self,options);
            end
        end
        
        function self = solve(self,config,specs,vars)
            % Compute plant state-space
            specs = rescale(specs,'constr+bound');
            [P,~,~,ch] = Solver.plant(config,specs,vars);
          
            % Weights for the norms in the objective
            alpha = zeros(1,length(specs.performance)); 
            alpha(1,1:specs.nobj) = cell2mat(cellfun(@(x) scale(x), specs.performance(1:specs.nobj), 'un', 0))';
            
            % Build the channel structure
            channels(length(ch.In)) = struct('performance',[],'In',[],'Out',[]);
            for j = 1:length(ch.In)
                if ismember(j,ch.Hinf)
                    channels(j).performance = inf;
                elseif ismember(j,ch.H2)
                    channels(j).performance = 2;
                else
                    error('cannot process performance type: contact the developers');
                end
                [channels(j).In,~] = find(ch.In{j});
                [~,channels(j).Out] = find(ch.Out{j});
            end
            % Solver input
            ncont = length(specs.ctrl_in);
            nmeas = length(specs.ctrl_out);
            P = balreal(std(P));
            [self.K,self.info] = LTI_unstructured_mix(P,nmeas,ncont,alpha,channels,specs.region,self.options);  
            % Save solver output
            self.performance = specs.performance;
            self.gamma = zeros(1,length(self.performance));
            self.mu = zeros(1,length(self.performance));
            %self.objective = self.info.objective;
            for k = 1:length(channels)
                if k <= specs.nobj
                    val = sqrt(self.info.gamma);
                else
                    val(k) = 1;
                end
                
                if channels(k).performance == Inf
                    self.gamma(k) = val(k);
                else
                    self.mu(k) = val(k);
                end
            end
                gamormu = self.gamma + self.mu;
                if specs.nobj > 0
                    obj = 1:specs.nobj;
                    self.performance(obj) = cellfun(@(x,y) {x*(1/y)}, self.performance(obj), num2cell(gamormu(1:specs.nobj)'));
                end
            end
    end
    
    methods (Static)        
        function cap = capabilities()
        % Returns the capabilities of \c mixedpoleMIMO. 
        % 
        % Return values: 
        %  cap: capabilities of \c mixedpoleMIMO @type struct
            cap.inout = 2;
            cap.norm = [2;Inf];
            cap.constraints = true;
            cap.unstable = false;
            cap.improper = false;
            cap.parametric = false;
            cap.fixedorder = false;
            cap.polereg = true;
        end
    end
end
