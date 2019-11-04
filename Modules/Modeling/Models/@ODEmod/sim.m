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

function [y,t,x] = sim(self,u,varargin)
    % Simulate the response of a system
    % Inputs
    % * self: the ODEmod
    % * u: function handle u(t) representing the input
    % * p: cell of parameter trajectories pj(t) or p(t) for single
    % parameter system
    % * t: simulation time. Either a time grid or the end time
    % * x0: initial state of the simulation
    % * key-value pairs: 
    %     - 'verbose': 0 or 1

    % event function for termination in steady state
    function [x,isterm,dir] = eventfun(t,x)
        dx = F(t,x);
        x = norm(dx)/norm(x) - 1e-3;
        isterm = 1;
        dir = 0;
    end
    % verbose function to print progress
    function status = outfun(t,y,flag)
        disp(t);
        status = 0;
    end

    % parse inputs
    % key-value pairs
    verbose = 0;
    isstr = find(cellfun(@ischar,varargin));
    if ~isempty(isstr)    
        keys = varargin(isstr);
        values = varargin(isstr+1);
        varargin([isstr,isstr+1]) = [];
        
        % set variables
        verbose = strcmp('verbose',keys);
        verbose = values{verbose};
    end
    
    iscel = cellfun(@iscell,varargin);
    issp = cellfun(@(x)isa(x,'function_handle'),varargin);
    isnum = cellfun(@isnumeric,varargin);
    assert(sum(isnum) <= 2,'Only 2 numerical inputs should be provided: simulation time and initial state');
    assert(sum(iscel) <= 1,'Only 1 cell should be provided: the scheduling parameters and their trajectory');
    assert(sum(issp) <= 1,'Only 1 function handle should be provided: the input trajectory');

    % inputs
    siz = size(u(0));
    assert(siz(1) == self.nu_, ['The number of input signals you are providing (' num2str(siz(1)) ') does not match the number of inputs of the system (' num2str(self.nu_) ').']);
    
    % parameters
    if ~any(iscel) || isempty(varargin{iscel})
        assert(~self.isparametric(), 'The system you are trying to simulate contains one or more scheduling parameter(s). Therefore, you have to define a trajectory for your parameters.');
        p = cell(0,2);
    else
        p = varargin{iscel};
        assert(all(cellfun(@(x) isa(x,'SchedulingParameter'), p(:,1))) && all(cellfun(@(x) isa(x,'function_handle'), p(:,2))), 'The scheduling parameter cell should contain SchedulingParameter objects in the first column and function handles in the second column.');
        assert(all(ismember(cellfun(@argument, self.parameters(), 'un', 0), cellfun(@argument, p(:,1),'un',0))), 'You didn''t define a trajectory for all scheduling parameters that are in the system you are trying to simulate.');
    end

    % time and initial state
    evt = false;
    switch(sum(isnum))
        case 0
            evt = true;
            t = [0, 100];
            x0 = self.x0;
        case 1
            t = varargin{isnum};
            x0 = self.x0;
        case 2
            [t,x0] = varargin{isnum};
            if x0~=self.x0; warning('You defined the initial states of your model earlier. I will ignore these for the simulation and use the ones you provided now.'); end
            assert(length(x0) == self.nx_, ['The length of your initial state vector (' num2str(length(x0)) ') does not match the number of states of your system (' num2str(self.nx_) ').']);
    end
    if length(t) == 1, t = [0,t]; end
    
    % Solve for the initial state vector in case some parts are not known yet  
    if any(isnan(x0))
        if self.Ts==0; u0 = u(t(1)); else; u0 = u(0); end
        
        disp('Solving for the internal states...');
        D_unknown = zeros(length(x0), sum(isnan(x0))); 
        D_unknown(sub2ind(size(D_unknown), find(isnan(x0))', 1:sum(isnan(x0)))) = 1;
        D_known = zeros(length(x0), sum(~isnan(x0))); 
        D_known(sub2ind(size(D_known), find(~isnan(x0))', 1:sum(~isnan(x0)))) = 1;
        nleq = @(x0_internal) self.g_conn(D_known*x0(~isnan(x0))+D_unknown*x0_internal, u0, cell2mat(cellfun(@(z)z(t), p(:,2), 'un', 0))) - x0_internal;
        x0_internal = fsolve(nleq, zeros(sum(isnan(x0)),1));
        x0(isnan(x0)) = x0_internal;
    end

    if self.Ts==0 
        % CONTINUOUS TIME %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Construct time dependent ODE
            if self.nparameters() > 0
                M = @(t,x) self.E(x,u(t),cell2mat(cellfun(@(z)z(t),p(:,2),'un',0)));
                F = @(t,x) self.f(x,u(t),cell2mat(cellfun(@(z)z(t),p(:,2),'un',0)));
            else
                M = @(t,x) self.E(x,u(t),[]);
                F = @(t,x) self.f(x,u(t),[]);
            end

            isODE = ~isempty(F(eps, x0));
            if isODE
                % set ODE options
                options = odeset('Mass',M,'RelTol',1e-4);
                if verbose, options = odeset(options,'OutputFcn',@outfun); end
                if evt
                    options.Events = @eventfun;
                end

                % solve the ODE
                disp('Solving the ODE or DAE...');
                [t,x] = ode23t(F,t,x0,options);
            %   [t,x] = ode15s(F,t,x0,options);
                if self.nparameters() > 0
                    y = cellfun(@(tau,xi)self.g(xi,u(tau),cell2mat(cellfun(@(x) x(tau), p(:,2), 'un',0))),num2cell(t),num2cell(x',1)','un',0);
                else
                    y = cellfun(@(tau,xi)self.g(xi,u(tau),0),num2cell(t),num2cell(x',1)','un',0);
                end
                y = transpose(cell2mat(transpose(y)));

            else
                % evaluate the output equation
                warning('You are simulating a static system, which is usually quite boring. Are you sure you didn''t forget anything?');
                if length(t)==2
                    t = linspace(t(1),t(2),100)';
                end
                y = transpose(cell2mat(cellfun(@(tau)self.g(zeros(0,1),u(tau),cell2mat(cellfun(@(z)z(tau),p(:,2),'un',0))),num2cell(t),'un',0)));
                x = zeros(length(t),0);
            end
    else
        % DISCRETE TIME %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        nsamp = ceil(diff(t)/self.Ts);
        x(1,:) = x0';
        for k = 1:nsamp
            y(k,:) = self.g(x(k,:)', u(k), cell2mat(cellfun(@(z)z(k),p(:,2),'un',0)))';
            t(k) = t(1)+k*self.Ts;
            if k==nsamp; break; end
            x(k+1,:) = self.f(x(k,:)',u(k),cell2mat(cellfun(@(z)z(k),p(:,2),'un',0)))';
        end
        
    end
end