%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build multivariate Lyapunov matrix and its time derivative
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUTS
%   basis           : TensorBasis
%   param           : cell with SchedulingParameters 
%   n               : dimension of Lyapunov matrix
%   opti            : to generate spline functions
%   dependency      : can create constant lyapunov matrices or matrices
%   dependent on the bounded rate of variations. It can be either
%   'constant' or 'brv.'
%   c_or_d          : to provide detaile to the function if it is for
%   continuous time LPV synthesis or discrete time. 'c' or 'd'.
%   varargin{1}     : number of degree of parameter
%   varargin{2}     : number of knots of parameter
%
% OUTPUTS
%   P    : [splines.Function] parameter-dependent Lyapunov matrix
%   dPdt : [splines.Function] time derivative of Lyapunov matrix


function [P,dPdt] = build_Lyapmat_new(basis,param,n,opti,dependency,c_or_d,varargin)
%
%   References :
%
%   P. Apkarian and R. J. Adams.
%   Advanced Gain-Scheduling Techniques for Uncertain Systems
%   IEEE Transactions on Control Systems Technology. 1998;6(1):21-32

check_optispline;
import splines.*;
np = length(param);
% construct Lyapunov matrix
switch c_or_d
    case 'c'
            j = 1;
            P_args = {};
            P_basis = {};
            for i = 1:np
                switch param{i}.bounded
                    case 0
                        switch dependency
                            case 'constant'
                                j = j + 1;
                            case 'brv'
                                P_args{j} = param{i}.tensor_basis.arguments;
                                P_basis{j} = basis.bases{i};
                                j = j + 1;
                            otherwise
                                error('Wrong selection of dependency')
                        end
                    case 'bounded'
                                P_args{j} = param{i}.tensor_basis.arguments;
                                P_basis{j} = basis.bases{i};
                                j = j + 1;
                    case 'unbounded'
                        j = j + 1;
                    otherwise
                end                
            end
            P = opti.Function(TensorBasis(P_basis,P_args),[n,n],'symmetric');
            
%                 for i = 1:length(param)
%                     switch 
%                 end
            % construct derivative
            switch dependency
                case 'constant'
                    dPdt = build_derivative(P);
                case 'brv'
                    dPdt = build_derivative(P,param);  
            end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lyapunov matrices for discrete time
    case 'd'
        switch dependency
            case 'constant'
                for i = 1:length(param)
                            P_basis{i} = basis.bases{i};
                            P_args{i}  = param{i}.tensor_basis.arguments;
                end
                     P = opti.Function(TensorBasis(P_basis,P_args),[n,n],'symmetric');
                     dPdt = P;
            case 'brv'
                error('Not yet implemented');
            otherwise
                error('Wrong selection of dependecy');
        end
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
    if (strcmp(param{i}.bounded,'bounded')) % if rate is not zero or unbounded, take derivative
        p = SchedulingParameter(strcat(param{i}.tensor_basis.arguments,'dot'),param{i}.rate,0);
        dPdt = dPdt + P.derivative(1,param{i}.tensor_basis.arguments)*p;
    end
end
end
end




