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


function [P,dPdt] = build_Lyapmat(basis,param,n,opti,dependency,c_or_d,varargin)
%
%   References :
%
%   P. Apkarian and R. J. Adams.
%   Advanced Gain-Scheduling Techniques for Uncertain Systems
%   IEEE Transactions on Control Systems Technology. 1998;6(1):21-32

check_optispline;
import splines.*;

% construct Lyapunov matrix
switch c_or_d
    case 'c'
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lyapunov matrices for discrete time
    case 'd'
            j = 1;
            P_args = {};
            P_basis = {};
        switch dependency
            case 'constant'
                for i = 1:length(param)
                        if  all(param{i}.rate == 0) % if rate is zero, add parameter
                            P_basis{j} = basis.bases{i};
                            bound = 0;
                            P_args{j}  = param{i}.tensor_basis.arguments;
                            j = j + 1;
                        else 
                            bound = 1;
                        end
                end
                     P = opti.Function(TensorBasis(P_basis,P_args),[n,n],'symmetric');
                     dPdt = P;
            case 'brv'

% This is the special case for discrete time LPV controller synthesis. This
% section allows to create Lyapunov matrices in discrete - time when the 
% scheduling parameters of the system are rate bounded. The main concept 
% of this section is to obtain a hexagonal mapping of the hyperrectangle 
% domain. It assumes that the desired knots of the scheduling parameters 
% are zero.

% References:
%
% J. De Caigny, J. F. Camino, R. C. L. F Oliveira, P. L. D. Peres and J.
% Swevers
% Gain-scheduled dynamic output feedback control for discrete-time LPV
% systems.
% International journal of Robust and Nonlinear Control (2012); 22:535-558

                        P = opti.Function(TensorBasis(P_basis,P_args),[n,n],'symmetric');
                        n_deg = varargin{1};
                        n_knots = varargin{2};
                        A = {};
                        dA = {};
                        Param_A = {};
                    for i = 1:length(param)
                        param_added{j} = param{i};
                        if ~all(param{i}.rate == 0) 
                        p = SchedulingParameter(strcat('d',param{i}.tensor_basis.arguments),param{i}.rate,0);
                        P_basis{j} = basis.bases{i};
                        dP_basis{j} = basis.bases{i};
                        % take for degree one and elevate it later after
                        % mapping
                            if ~(n_deg(i)==1)
                                for k = n_deg(i):-1:2
                                P_basis{i} = P_basis{i}.derivative;
                                dP_basis{i} = dP_basis{i}.derivative;
                                end
                            end
                            if ~any(isinf(param{i}.rate)) % if rate is not unbounded, add parameter
                            P_args{j}  = param{i}.tensor_basis.arguments;
                            l{i} = param{i}.basis.domain.min;
                            u{i} = param{i}.basis.domain.max;
                            rate{i} = param{i}.rate;
                            dl{i} = rate{i}(1);
                            du{i} = rate{i}(2);
                            del = 0.001;   
                                 P_basis{j} = P_basis{j}.insert_knots([l{i}-dl{i}+del,u{i}-du{i}-del]);
                                 dP_basis{j} = dP_basis{j}.insert_knots([l{i}-dl{i},u{i}-du{i}]);
                                 Param_A = [Param_A, {param_added{j}.tensor_basis.arguments}];
                                 j = j + 1;
                                 P_args{j} = p.tensor_basis.arguments;
                                 P_basis{j} = BSplineBasis([p.basis.domain.min*ones(1,1);linspace(p.basis.domain.min,p.basis.domain.max,0+2)';p.basis.domain.max*ones(1,1)],1);
                                 dP_basis{j} = BSplineBasis([p.basis.domain.min*ones(1,1);linspace(p.basis.domain.min,p.basis.domain.max,0+2)';p.basis.domain.max*ones(1,1)],1);
                                 param_added{j} = p;
                            else
                                j = i + 1;
                            end
                %rectangle to hexagon mapping
 
                                H{i} = [l{i},    l{i},          l{i}-dl{i}+del,     u{i}-du{i}-del,     u{i},    u{i},        u{i}-du{i}-del,        l{i}-dl{i}+del;
                                        0,       du{i}+del,     du{i}+del,          du{i}+del,          0,      dl{i}-del,    dl{i}-del,             dl{i}-del];

                                Hb{i} = [l{i},   l{i},      l{i},          l{i},           l{i}-dl{i},      l{i}-dl{i}+del,        u{i}-du{i},     u{i}-du{i}-del,  u{i},       u{i},       u{i},     u{i},        u{i}-du{i},      u{i}-du{i}-del,        l{i}-dl{i},         l{i}-dl{i}+del;
                                        0,       0,         du{i},          du{i}+del       du{i},           du{i}+del,             du{i},           du{i}+del,     0,          0,      dl{i},     dl{i}-del,   dl{i},           dl{i}-del,             dl{i},              dl{i}-del];

                                Cx{i} = [H{i}(1,1), H{i}(1,8), H{i}(1,7), H{i}(1,6);
                                        H{i}(1,2), H{i}(1,3), H{i}(1,4), H{i}(1,5)]';

                                Cy{i} = [H{i}(2,1), H{i}(2,8), H{i}(2,7), H{i}(2,6);
                                        H{i}(2,2), H{i}(2,3), H{i}(2,4), H{i}(2,5)]';

                                Cxb{i} = [Hb{i}(1,1), Hb{i}(1,2), Hb{i}(1,15), Hb{i}(1,16), Hb{i}(1,13), Hb{i}(1,14), Hb{i}(1,11), Hb{i}(1,12);
                                          Hb{i}(1,3), Hb{i}(1,4), Hb{i}(1,5), Hb{i}(1,6), Hb{i}(1,7), Hb{i}(1,8), Hb{i}(1,9), Hb{i}(1,10)]';

                                Cyb{i} = [Hb{i}(2,1), Hb{i}(2,2), Hb{i}(2,15), Hb{i}(2,16), Hb{i}(2,13), Hb{i}(2,14), Hb{i}(2,11), Hb{i}(2,12);
                                         Hb{i}(2,3), Hb{i}(2,4), Hb{i}(2,5), Hb{i}(2,6), Hb{i}(2,7), Hb{i}(2,8), Hb{i}(2,9), Hb{i}(2,10)]';
                                 a{i} = Function(TensorBasis({P_basis{j-1:j}}),Cx{i});
                                 da{i} = Function(TensorBasis({dP_basis{j-1:j}}),Cy{i});
                                 figure;
                                  n_l = 101;
                                x= linspace(l{i},u{i},n_l)';
                                y= linspace(dl{i}-del,du{i}+del,n_l)';
                                af = a{i}.grid_eval({x,y(:,1)});
                                daf = da{i}.grid_eval({x,y(:,1)});
                                plot(af(:), daf(:), '.'), hold all
                                xlabel('p(k)');
                                ylabel('dp');
                                title('rectangle to hexagonal mapping of parameters');
                        A = [A,{a{i}.basis(1),a{i}.basis(2)}];
                        dA = [dA, {a{i}.basis(1)+da{i}.basis(1),a{i}.basis(2)+da{i}.basis(2)}];
                        Param_A = [Param_A, {param_added{j}.tensor_basis.arguments}];
                        else
                            dPdt = P;
                        end
                        j = j + 1;
                    end
                    P = opti.Function(TensorBasis(A,Param_A),[n,n],'symmetric');
                        for l = 1:length(param)
                            if ~(all(param{l}.rate==0))
                                dPdt = opti.Function (TensorBasis(dA,Param_A),[n,n],'symmetric');                                
                            else
                                dPdt = P;  
                            end
                        end
                        for l = 1:length(param)
                            if ~(all(param{l}.rate==0))
                                P = (a{l}^(n_deg(l)-1))*P;
                                dPdt = ((a{l}+da{l})^(n_deg(l)-1))*dPdt;
                            end
                        end
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
    if ~(all(param{i}.rate == 0)) && ~any(isinf(param{i}.rate)) % if rate is not zero or unbounded, take derivative
        p = SchedulingParameter(strcat(param{i}.tensor_basis.arguments,'dot'),param{i}.rate,0);
        dPdt = dPdt + P.derivative(1,param{i}.tensor_basis.arguments)*p;
    end
end
end
end




