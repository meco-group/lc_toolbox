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

%% TIME DOMAIN SIMULATIONS OF IMPROPER LPV SYSTEMS 
%  (for integration in LCToolbox)
%
%  This function tackles time domain simulations for the most general
%  system type available in LCToolbox, being improper LPV systems.
% 
%  IMPORTANT NOTE: If you know in advance that your system has a proper
%  representation with E = I for the whole parameter space, this solver is 
%  not very efficient (the calculations of the transformations are then 
%  obsolete). The derivatives of the 'improper states' are assumed to be
%  zero initially.
%
%  Laurens Jacobs, MECO Research Team, KU Leuven, May 2017

function varargout = lsimlpv(varargin)

import splines.*
%  ___________________
% | 1. INPUT PARSING  |
% |___________________|

% 1.1  Input arguments list

switch nargin
    case 4
        model = varargin{1};
        u = varargin{2};
        p = varargin{3};
        t = varargin{4};
    case 5
        model = varargin{1};
        u = varargin{2};
        p = varargin{3};
        t = varargin{4};
        x0 = varargin{5};
    otherwise
        error('Wrong number of input arguments.')
end

% 1.2 Tons of input checks and reordering of the parameters

assert(isa(model,'AbstractLPVmod'),'Your model is not an LPV model.');
assert(isnumeric(t) && range(diff(t)) <= 1e-10 && all(diff(t) > 0),'Your time vector should be equally spaced and monotonically increasing.');

if isnumeric(u)
    assert(size(u,1) == size(t,1),'The length of the input vectors and time vector should be equal.');
    uu = u;
    u = @(k) interp1(t,uu,k)';
    iter = length(t)-1;
elseif isa(u,'function_handle')
    iter = length(t)-1;
else
    error('Invalid type for the inputs: they can only be defined by a matrix or by a function handle.');
end

if all(cellfun(@isnumeric,p(:,2)))
    assert(all(cellfun(@(x) size(x,1) == size(t,1),p(:,2))),'The length of the parameter vectors and time vector should be equal.');
    pp = cell2mat(cellfun(@transpose,p(:,2),'un',0));
    for i = 1:size(p,1)
        p{i,2} = @(k) interp1(t,pp(i,:)',k)';
    end
    iter = length(t)-1;
elseif all(cellfun(@(x) isa(x,'function_handle'), p(:,2)))
else
    error('Invalid type for the parameters: they can only be defined by vectors or by function handles.');
end

assert(size(p,1) >= size(model.parameters(),2),'The cell with your parameters and their values or function handles does not contain all the parameters of the model.');
ordened_params = cell(size(model.parameters(),2),1);
for i = 1:size(model.parameters(),2)
    chk = false;
    for j = 1:size(p,1)
        chk = chk || isequal(p{j,1},model.parameters{1,i});
        if chk; ordened_params{i} = p{j,2}; break; end;
    end
    assert(chk,'The cell with your parameters and their values or function handles does not contain all the parameters of the model.');
end

if exist('x0','var')
    assert(all(size(x0) == [model.nx(), 1]),'The length of x0 is different from the number of states your system has.');
else
    if rank(model.E) ~= length(model.E)
        warning('Not providing initial values for the improper states consistent with the system dynamics sometimes leads to unexpected results.');
    end
    x0 = zeros(model.nx(),1);
end

%  __________________________________
% | 2. SOLVE THE IMPROPER ODE (RK4)  |
% |__________________________________|

% initialize
X = zeros(model.nx(),length(t));
Y = zeros(model.nout(),length(t));

% time step
h = diff(t(1:2));

% first time instance
i = 1;
[A4,B4,C4,D4,E4] = evalSysMat(model,p,t(i));
[A4p,B4p,~,~,~] = proper_formulation(A4,B4,C4,D4,E4);
xp = x0(1:rank(E4));
Y(:,i) = C4*x0 + D4*u(t(i));

for i = 2:iter
    
    % evaluate system matrices at time instances needed for RK4
    [A23,B23,C23,D23,E23] = evalSysMat(model,p,t(i)+h/2);
    [A4,B4,C4,D4,E4] = evalSysMat(model,p,t(i)+h);
    
    % evaluate the input at the time instances needed
    u1 = u(t(i-1));
    u23 = u(t(i-1)+h/2);
    u4 = u(t(i));

    % get rid of the improper part if present
    A1p = A4p; B1p = B4p;
    [A23p,B23p,~,~,~] = proper_formulation(A23,B23,C23,D23,E23);
    [A4p,B4p,F4,Tp2ip_x4,Tp2ip_u4] = proper_formulation(A4,B4,C4,D4,E4);

    % RK4 increments
        % increment 1
        x1 = xp;
        k1p = A1p*x1(1:length(A1p)) + B1p(1:length(A1p),:)*u1;
        
        % increment 2
        x2 = xp+h/2*k1p; 
        k2p = A23p*x2(1:length(A23p)) + B23p(1:length(A23p),:)*u23;
        
        % increment 3
        x3 = xp+h/2*k2p;
        k3p = A23p*x3(1:length(A23p)) + B23p(1:length(A23p),:)*u23;  
        
        % increment 4
        x4 = xp+h*k3p; 
        k4p = A4p*x4(1:length(A4p)) + B4p(1:length(A4p),:)*u4;
     
    % evaluate new proper 'states' (these are not the original states!)
    xp = xp + h/6*(k1p+2*k2p+2*k3p+k4p);
    
    % evaluate new improper states and store in the global state vector
    xip = Tp2ip_x4*xp + Tp2ip_u4*u4;
    
    % evaluate actual proper states
    Xp = xp - F4*xip;
    
    % evaluate output
    X(:,i) = [Xp ; xip];
    Y(:,i) = C4*X(:,i) + D4*u4;
    
end
    
    if i ~= length(t)
        Y(:,i+1) = Y(:,i);  % copy last value in case of numeric inputs
        X(:,i+1) = X(:,i);  % copy last value in case of numeric inputs
    end
    
    if nargout == 0 % graphical output
        figure;
        for i = 1:size(Y,1)
            h(i) = subplot(size(Y,1),1,i);
            plot(t,Y(i,:));
            ylabel(['To: output ' num2str(i)]);
            axis tight; 
            if i ~= size(Y,1); set(gca, 'XTickLabel', []); end;
        end
        xlabel('Time');
        h(i) = subplot(size(Y,1),1,1);
        title('LPV Simulation Results');
    else % numerical output
        varargout = {Y,t,X};
    end;
    
end

function [A,B,C,D,E] = evalSysMat(model,p,t)
    if isa(model,'AbstractDSSmod')
        [A,B,C,D,E] = getdssdata(model.evalme(cellfun(@feval,p(:,2),num2cell(repmat(t,length(p(:,2)),1))),model.parameters()'));
    else % isa(model,'AbstractLFTmod')
        [A,B,C,D,E] = getdssdata(model.evalme(cellfun(@feval,p(:,2),num2cell(repmat(t,length(p(:,2)),1))),model.parameters()'));
    end
end

function [Ap,Bp,F,Tp2ip_x,Tp2ip_u] = proper_formulation(A,B,C,D,E)
    
    % diagonalize E
    p = rank(E);
    [sys,~] = rrefE(dss(A,B,C,D,E));
    [A,B,~,~,E] = dssdata(sys);
    
    % separate the proper part of the dynamics from the improper part
    Ap = A(1:p,1:p); Bp = B(1:p,:);
    F = E(1:p,p+1:end);
    
    % substitution trick (working back towards the original dynamics)
    sub = rref([A(p+1:end,p+1:end) A(p+1:end,1:p) B(p+1:end,:)]);
    subA = -sub(:,end-p-size(B,2)+1:end-size(B,2));
    subB = -sub(:,end-size(B,2)+1:end);
    for i = 1:p
        for j = p+1:size(A,2)
            Ap(i,:) = Ap(i,:) + A(i,j)*subA;
            Bp(i,:) = Bp(i,:) + A(i,j)*subB;
        end
    end
    
    % provide the backtransformation matrices
    % from proper states to proper + improper states:
    % xip = Tp2ip_x * xp + Tp2ip_u * u
    Tp2ip_x = subA;
    Tp2ip_u = subB;

end