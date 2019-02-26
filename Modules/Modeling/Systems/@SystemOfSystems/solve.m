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

function [config,K,solver] = solve(config,obj,constr,varargin)
%SOLVE Summary of this function goes here
%   Detailed explanation goes here

% 1. Parse the inputs for variables and parameters
vars = {};
modelnames = {};
options = [];

for k = 1:length(varargin)
    if isa(varargin{k},'SystemOfModels')
        vars = [vars,varargin(k)];
    elseif ischar(varargin{k})
        modelnames = [modelnames;varargin(k)];
    elseif isstruct(varargin{k})
        options = varargin{k};
    end
end

if isempty(vars)
    vars = get_empty(config);
end

% 2. Construct a controller design problem
specs = ControllerDesign();
% objectives
if ~isa(obj, 'cell') && ~isempty(obj); obj = {obj}; elseif isempty(obj); obj = {}; end
assert(all(cellfun(@(x) isa(x,'Norm'), obj)), 'Objectives can only be of type Norm.');
assert(size(obj,2)<=1, 'Please stack objectives only vertically.');
specs = specs.addobjective(obj);
% constraints
if ~isa(constr, 'cell') && ~isempty(constr); constr = {constr}; elseif isempty(constr); constr = {}; end
assert(all(cellfun(@(x) isa(x,'NormConstraint') || isa(x,'Order') || isa(x,'Region'), constr)),'Constraints can only be of type NormConstraint, Order or Region.');
assert(size(constr,2)<=1, 'Please stack constraints only vertically.');
orders = cellfun(@(x)isa(x,'Order'),constr);
regions = cellfun(@(x)isa(x,'Region'),constr);
assert(sum(orders)<=1,'Multiple order constraints are not allowed.');
if any(orders), specs.order = order(constr{orders}); constr(orders) = []; end
if any(regions), specs.region = constr(regions); constr(regions) = []; end
specs = specs.addconstraint(constr);
% control in and outputs
ctrl_in_ = cellfun(@(x)x.out,vars,'UniformOutput',false);
specs.ctrl_in = vertcat(ctrl_in_{:});
ctrl_out_ = cellfun(@(x)x.in,vars,'UniformOutput',false);
specs.ctrl_out = vertcat(ctrl_out_{:});

% 4. Choose an appropriate solver and solve
[config,r] = subsystem(config,modelnames{:});
solver = Solver.select(config,specs,vars,options);
solver = solver.solve(config,specs,vars);
restore(r);

% 5. Add the solution to the variable
assert(length(vars) == 1,'Solving multiple controllers at once not yet implemented');
K = solver.K;
if ~isfield(options,'controller_name')
    K.name = ['controller' num2str(vars{1}.numod())];
else    
    K.name = options.controller_name;
end
vars{1}.add(K);
solver.K.name = K.name;

% 6. Add designed closed loop and specs to solver object
[~,r] = vars{1}.empty();
vars{1}.add(solver.K);
solver.H = config.model();
vars{1}.empty();
restore(r);

end

