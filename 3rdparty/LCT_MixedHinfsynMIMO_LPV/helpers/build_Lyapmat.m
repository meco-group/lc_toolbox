%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build multivariate Lyapunov matrix and its time derivative
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUTS
%   basis : TensorBasis
%   param : cell with SchedulingParameters 
%   n     : dimension of Lyapunov matrix
%   opti  : OptiSplineYalmip function
%   dependency : If the dependency to be selected for Lyapunov matrices is
%   constant or bounded rate of variations
%
% OUTPUTS
%   P    : [splines.Function] parameter-dependent Lyapunov matrix
%   dPdt : [splines.Function] time derivative of Lyapunov matrix
%
%   References :
%
%   P. Apkarian and R. J. Adams.
%   Advanced Gain-Scheduling Techniques for Uncertain Systems
%   IEEE Transactions on Control Systems Technology. 1998;6(1):21-32
   
function [P,dPdt] = build_Lyapmat(basis,param,n,opti,dependency)

import splines.*

% construct Lyapunov matrix
j = 1;
P_args = {};
P_basis = {};
switch dependency
    case 'constant'
        for i = 1:length(param)
            if (all(param{i}.rate == 0) && ~any(isinf(param{i}.rate))) % if rate is zero, add parameter
                P_args{j}  = param{i}.tensor_basis.arguments;
                P_basis{j} = basis.bases{i};
                j = j + 1;
            else 
            end
        end
    case 'brv'
        for i = 1:length(param)
            if ~any(isinf(param{i}.rate)) % if rate is not unbounded, add parameter
                P_args{j}  = param{i}.tensor_basis.arguments;
                P_basis{j} = basis.bases{i};
                j = j + 1;
            end
        end
end
P = opti.Function(TensorBasis(P_basis,P_args),[n,n],'symmetric');

% construct derivative
switch dependency
    case 'constant'
        dPdt = build_derivative(P);
    case 'brv'
        dPdt = build_derivative(P,param);  
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ dPdt ] = build_derivative(P,varargin)
% This function computes the derivative of any B-spline function and
% returns zeros if the variable is constant and returns rate dependent
% function in case of varying variable.
n = size(P,1);
dPdt = zeros(n);
if (nargin == 2)
    param = varargin{1};
for i = 1:length(param)
    if ~(all(param{i}.rate == 0)) && ~any(isinf(param{i}.rate)) % if rate is not zero or unbounded, take derivative
        p = SchedulingParameter(strcat(param{i}.tensor_basis.arguments,'dot'),param{i}.rate,0);
        dPdt = dPdt + P.derivative(1,param{i}.tensor_basis.arguments)*p;
    end
end
end
end




