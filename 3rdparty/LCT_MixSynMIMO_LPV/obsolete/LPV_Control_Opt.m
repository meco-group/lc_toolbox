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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to specify options for the LPV control design
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function self = LPV_Control_Opt(varargin)

required = {};
optional = {'spec','objective','var_deg','var_knots','relax_deg','relax_mp','tolerance','verbose','scaling_obj','solver','controller_dependency'};

calls = varargin(1:2:end);

assert(mod(nargin,2) == 0,'Odd number of arguments.');
assert(iscellstr(calls),'The names of your name/value pairs should be strings.');
assert(~any(~ismember(required,calls)),['These arguments are required: ' strjoin(required)]);

for i = 1:nargin/2
    if strcmp(varargin{2*i-1},'spec'); assert(~exist('spec','var'),'You cannot define the performance specification twice.'); spec = varargin{2*i}; end;
    if strcmp(varargin{2*i-1},'objective'); assert(~exist('objective','var'),'You cannot define the objective function twice.'); objective = varargin{2*i}; end;
    if strcmp(varargin{2*i-1},'var_deg'); assert(~exist('var_deg','var'),'You cannot define the polynomial degree of the optimization variables twice.'); var_deg = varargin{2*i}; end;
    if strcmp(varargin{2*i-1},'var_knots'); assert(~exist('var_knots','var'),'You cannot define the number of internal knots of the optimization variables twice.'); var_knots = varargin{2*i}; end;
    if strcmp(varargin{2*i-1},'relax_deg'); assert(~exist('relax_deg','var'),'You cannot define the relaxation degree twice.'); relax_deg = varargin{2*i}; end;
    if strcmp(varargin{2*i-1},'relax_mp'); assert(~exist('relax_mp','var'),'You cannot define the number of midpoint refinements twice.'); relax_mp = varargin{2*i}; end;
    if strcmp(varargin{2*i-1},'tolerance'); assert(~exist('tolerance','var'),'You cannot define the tolerance of the primal residual twice.'); tolerance = varargin{2*i}; end;
    if strcmp(varargin{2*i-1},'verbose'); assert(~exist('verbose','var'),'You cannot define verbose twice.'); verbose = varargin{2*i}; end;
    if strcmp(varargin{2*i-1},'scaling_obj'); assert(~exist('scaling_obj','var'),'You cannot define a scaling for the objective function twice.'); scaling_obj = varargin{2*i}; end;
    if strcmp(varargin{2*i-1},'solver'); assert(~exist('solver','var'),'You cannot define an SDP solver twice.'); solver = varargin{2*i}; end;
    if strcmp(varargin{2*i-1},'controller_dependency'); assert(~exist('controller_dependency','var'),'You cannot define controller dependency twice.'); controller_dependency = varargin{2*i}; end;
    if ~ismember(varargin{2*i-1},[required, optional]); warning(['I could not identify the option ''' varargin{2*i-1} ''' and I am going to ignore it.']); end;
end

if ~exist('spec','var'); spec = inf; end;
switch spec
    case 0
        objective = casadi.MX(0);
        tolerance = 0;
    case {2,inf}
        if ~exist('objective','var'); objective = 'wc'; end;
        if ~exist('tolerance','var'); tolerance = 1e-6; end;
end
if ~exist('var_deg','var'); var_deg = 1; end;
if ~exist('var_knots','var'); var_knots = 0; end;
if ~exist('relax_deg','var'); relax_deg = 0; end;
if ~exist('relax_mp','var'); relax_mp = 0; end;
if ~exist('verbose','var'); verbose = 0; end;
if ~exist('scaling_obj','var'); scaling_obj = 1; end;
if ~exist('solver','var'); solver = 'mosek'; end;
if ~exist('controller_dependency','var'); controller_dependency = 'a'; end;

self = LPV_Control_Opt1(spec,objective,var_deg,var_knots,relax_deg,relax_mp,tolerance,verbose,scaling_obj,solver,controller_dependency);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
function disp(self)

disp([' Performance specification: ' num2str(self.spec)]);
if ischar(self.objective)
    disp([' Objective function: ' num2str(self.objective)]);
else
    disp([' Objective function: ' num2str(0)]);
end
disp([' Degree optimization variables: ' num2str(self.var_deg)]);
disp([' Number of knots optimization variables: ' num2str(self.var_knots)]);
disp([' Degree LMI relaxation: ' num2str(self.relax_deg)]);
disp([' Midpoint refinements LMI relaxation: ' num2str(self.relax_mp)]);
disp([' Tolerance primal residual SDP solver: ' num2str(self.tolerance)]);
disp([' Verbose: ' num2str(self.verbose)]);
disp([' Scaling of objective function: ' num2str(self.scaling_obj)]);
disp([' Solver: ' num2str(self.solver)]);
disp([' Controller parameter dependency: ' num2str(self.controller_dependency)]);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
function self = LPV_Control_Opt1(spec,objective,var_deg,var_knots,relax_deg,relax_mp,tolerance,verbose,scaling_obj,solver,controller_dependency)

p = inputParser;
addRequired(p,'spec',@(x) assert(isnumeric(spec),'Performance specification should be 0, 2 or inf, corresponding to the design of a stabilizing, H-2 or H-infinity controller).'));
addRequired(p,'objective',@(x) assert(ischar(objective) || isa(objective,'casadi.MX'),'The objective should be a string (wc,L1).'));
addRequired(p,'var_deg',@(x) assert(isvector(var_deg),'The degree of the optimization variables should be a scalar or vector whose length corresponds to the number of parameters.'));
addRequired(p,'var_knots',@(x) assert(isvector(var_knots),'The number of knots of the optimization variables should be a scalar or vector whose length corresponds to the number of parameters.'));
addRequired(p,'relax_deg',@(x) assert(isvector(relax_deg),'The relaxation degree should be a scalar or vector whose length corresponds to the number of parameters.'));
addRequired(p,'relax_mp',@(x) assert(isvector(relax_mp),'The number of midpoint refinements should be a scalar or vector whose length corresponds to the number of parameters.'));
addRequired(p,'tolerance',@(x) assert(isscalar(tolerance),'The tolerance of the primal residual should be a scalar.'));
addRequired(p,'verbose',@(x) assert(isscalar(verbose),'Verbose should be a scalar.'));
addRequired(p,'scaling_obj',@(x) assert(isscalar(scaling_obj),'The objective scaling should be a scalar.'));
addRequired(p,'solver',@(x) assert(ischar(solver),'SDP solver should be a string (mosek,sdpt3,...)'));
addRequired(p,'controller_dependency',@(x) assert(ischar(controller_dependency),'Controller dependency should be a string (a,ada)'));

parse(p,spec,objective,var_deg,var_knots,relax_deg,relax_mp,tolerance,verbose,scaling_obj,solver,controller_dependency);

% assign properties
self.spec = p.Results.spec;
self.objective = p.Results.objective;
self.var_deg = p.Results.var_deg;
self.var_knots = p.Results.var_knots;
self.relax_deg = p.Results.relax_deg;
self.relax_mp = p.Results.relax_mp;
self.tolerance = p.Results.tolerance;
self.verbose = p.Results.verbose;
self.scaling_obj = p.Results.scaling_obj;
self.solver = p.Results.solver;
self.controller_dependency = p.Results.controller_dependency;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   