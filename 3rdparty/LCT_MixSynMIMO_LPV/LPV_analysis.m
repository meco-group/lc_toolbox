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
% ANALYSIS FOR LPV SYSTEMS USING BSPLINES 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analysis for the MIMO LPV system
%   dx =  A(a)x +  B(a)w
%   z  = C(a)x + D(a)w
% we consider the design of an unstructured output-feedback controller
%   dxc = Ac(a,da)xc + Bc(a)y
%     u =    Cc(a)xc + Dc(a)y
% Hence, the resulting controller is of the same order as the generalized
% plant, and depends on the parameter and (optionally) its rate of variation. 
%
% -> continuous time
% -> state-space matrices have a B-spline dependency on multivariate
% parameter, where each scalar parameter:
%    1) is constant 
%    2) has a bounded rate of variation
%    3) has an unbounded rate of variation
%
% Performance objective:
% -> stability, H-2 norm or H-infinity norm
% -> worst-case or L1 norm optimization 
%
% LMI relaxations:
% -> midpoint refinements
% -> degree elevation
%
% INPUT ARGUMENTS
% sys        : instance of ... class
% ny         : dimension measured output
% nu         : dimension control input
% param      : cell with SchedulingParameters
% opti       : something Casadi related...
% opts       : instance of LPV_Control_Opt class
%
% OUTPUT ARGUMENTS
% solution : struct with primal solution and details
%
% Uses the cpp_spline toolbox available online: 
%   https://gitlab.mech.kuleuven.be/meco-software/cpp_splines
%
% Taranjitsingh Singh, Gijs Hilhorst, Div. PMA, Dept. Mechanical Engineering, KU Leuven, 2017
%
%
% References:
%   
%   P. Apkarian and R. J. Adams.
%   Advanced Gain-Scheduling Techniques for Uncertain Systems
%   IEEE Transactions on Control Systems Technology. 1998;6(1):21-32
%
%   C.W. Scherer, P. Gahinet and M. Chilali.
%   Multiobjective Output-Feedback Control via LMI Optimization.
%   IEEE Transactions on Automatic Control. 1997;42(7):896-911   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Primal,Dual] = LPV_analysis(sys,param,varargin)

default = struct('spec',inf,'objective','L1','var_deg',2,...
    'var_knots',0,'relax_deg',0,'relax_mp',0,'tolerance',1e-6,...
    'verbose',0,'scaling_obj',1,'solver','mosek',...
    'controller_dependency','a');

if nargin == 3
    opts = mergestruct(varargin{1},default);
else
    opts = default;
end

% solve primal problem
t = cputime;
Primal = LPV_analysis_primal(sys,param,opts);
Primal.tc = cputime - t;
if isempty(Primal)
    Dual = []; return;
end
Dual = [];
% % solve dual problem
%t = cputime;
% Dual = LPV_analysis_dual(sys,ny,nu,param,opts);
% Dual.tc = cputime - t;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sol_info = LPV_analysis_primal(sys,param,opts)

check_optispline;
import splines.*;

t_start = cputime; 
options = sdpsettings('solver',opts.solver,'verbose',opts.verbose);
opti = OptiSplineYalmip();
gen_sys = extract_generalized_plant(sys,0,0);
%LPV System
A = gen_sys.A; nx = gen_sys.nx;
B = sys.B; C = sys.C; D = sys.D; nw = gen_sys.nw; nz = gen_sys.nz;

N = length(param); % take number of parameters

% determine degree and knots for each basis
if     length(opts.var_deg) == 1; deg = opts.var_deg*ones(1,N);
elseif length(opts.var_deg) == N; deg = opts.var_deg;
else   error('opts.var_deg has incorrect dimension'); end
if     length(opts.var_knots) == 1; nknots = opts.var_knots*ones(1,N);
elseif length(opts.var_knots) == N; nknots = opts.var_knots;
else   error('opts.var_knots has incorrect dimension'); end

for i = 1:N
    args{i} = param{i}.tensor_basis.arguments;
    Bt{i} = BSplineBasis([param{i}.basis.domain.min*ones(deg(i),1);linspace(param{i}.basis.domain.min,param{i}.basis.domain.max,nknots(i)+2)';param{i}.basis.domain.max*ones(deg(i),1)],deg(i));
end
TB = TensorBasis(Bt,args);
switch sys.Ts
    case {0,{}}
        [P,dPdt] = build_Lyapmat(TB,param,nx,opti,'brv','c');
        Phi = [dPdt*eye(nx), P; P, zeros(nx)];
    otherwise
        [P,dPdt] = build_Lyapmat(TB,param,nx,opti,'brv','d',deg,nknots);
        Phi = [-P, zeros(nx); zeros(nx), dPdt*eye(nx)];
end
    switch opts.spec
       case inf
            switch opts.objective
                case 'wc'; gam2 = opti.variable(1);
                case 'L1'; gam2 = opti.Function(TB,[1,1],'symmetric');
                otherwise; error('PRIMAL SDP: invalid performance objective');
            end
        case 2
            W = opti.Function(TB,[nz,nz],'symmetric');
            switch opts.objective
                case 'wc'; gam2 = opti.variable(1);
                case 'L1'; gam2 = opti.Function(TB,[1,1],'symmetric');
                otherwise; error('PRIMAL SDP: invalid performance objective');
            end
        otherwise
            error('wrong selection of specification');
    end

    switch opts.spec
        case inf
%             T11 = A'*P + P*A + dPdt;
%             T21 = B'*P;
%             T22 = -gam2*eye(nw);
%             T31 = C;
%             T32 = D;
%             T33 = -gam2*eye(nz);
% 
%             Term = [T11, T21', T31';...
%                     T21, T22 , T32';...
%                     T31, T32 , T33];

        inner_factor = [Phi, zeros(2*nx,nw+nz); zeros(nw+nz,2*nx), [-gam2*eye(nw), zeros(nw,nz); zeros(nz,nw), eye(nz)]];
        outer_factor = [eye(nx)     , zeros(nx,nw); 
                        A           , B           ; 
                        zeros(nw,nx), eye(nw)     ;
                        C           , D          ];
        Term = outer_factor'*inner_factor*outer_factor;
                
        case 2
        inner_factor = [Phi, zeros(2*nx,nw); zeros(nw,2*nx), -eye(nw)];
        outer_factor = [eye(nx)     , zeros(nx,nw); 
                        A           , B           ; 
                        zeros(nw,nx), eye(nw)    ];
        Term1 = outer_factor'*inner_factor*outer_factor;
        Term2 = [W, C; C', P];
        Term3 = trace(W) - gam2;    
        Term = blkdiag(Term1,-Term2,Term3);
            
        otherwise
            error('wrong selection of specification');
    end    

% LMI relaxations

% degree elevations
if any(opts.relax_deg)
    if     length(opts.relax_deg) == 1; relax_deg = opts.relax_deg*ones(1,N);
    elseif length(opts.relax_deg) == N; relax_deg = opts.relax_deg;
    else   error('opts.relax_deg has incorrect dimension'); end
    Q    = Q.degree_elevation(relax_deg,args);
    Term = Term.degree_elevation(relax_deg,args);
end

% midpoint refinements
if any(opts.relax_mp)
    if     length(opts.relax_mp) == 1; relax_mp = opts.relax_mp*ones(1,N);
    elseif length(opts.relax_mp) == N; relax_mp = opts.relax_mp;
    else   error('opts.relax_mp has incorrect dimension'); end
    Q    = Q.midpoint_refinement(relax_mp,args);
    Term = Term.midpoint_refinement(relax_mp,args);
end

    
comp_time_parsing = cputime - t_start;
    
% solve SDP
constraints = {Term <= 0,P >= 0};
switch opts.spec
    case 0
        objective = casadi.MX(0);
    case {2,inf}
        switch opts.objective
            case 'wc'; objective = gam2;
            case 'L1'; objective = gam2.integral;
        end
        objective = objective*opts.scaling_obj;
end
opti.minimize(objective)
opti.subject_to(constraints)
opti.solver('yalmip',struct('yalmip_options',options));
sol = opti.solve();
comp_time_SDP = cputime - t_start - comp_time_parsing;

% extract solution info
sol_info.coeffs_P    = P.coeff.dimension;
sol_info.coeffs_Term = Term.coeff.dimension;
[pr,~]               = checkset(opti.yalmip_constraints); % primal residual
if pr > -opts.tolerance
    sol_info.feasible  = 1;
    sol_info.objective = sqrt(sol.value(objective));
    switch opts.spec
        case {2,inf}
            sol_info.gam2 = sol.value(gam2);
    end
    sol_info.P      = sol.value(P);
    sol_info.Term   = sol.value(Term);
else
    sol_info.feas      = 0;
    sol_info.objective = []; 
    disp('PRIMAL SDP: infeasible problem'); return;
end
sol_info.comp_time_parsing = comp_time_parsing;
sol_info.comp_time_SDP     = comp_time_SDP;
sol_info.comp_time_total   = cputime - t_start;

% check solution
for j = 1:prod(sol_info.P.coeff.dimension)
    min_eig(j) = min(eig(sol_info.P.coeff.data(j))); 
end
sol_info.check_pos_P = min(min_eig);
for j = 1:prod(sol_info.Term.coeff.dimension)
    max_eig(j) = max(eig(sol_info.Term.coeff.data(j))); 
end
sol_info.check_neg_Term = max(max_eig);    

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function sys = generalized_plant(sys,nu,ny)         
            if nargin == 1; sys.nu = 0; sys.ny = 0; end
            if nargin == 2; sys.ny = 0; end
            sys.nx = size(sys.A,1);
            sys.nu = nu;
            sys.nw = size(sys.B,2)-nu;
            sys.ny = ny;
            sys.nz = size(sys.C,1)-ny;
            if sys.nu ~= 0 && sys.ny == 0     % state feedback case
                sys.Bw = temporary_fix_subsref(sys.b,1:sys.nx,1:sys.nw); sys.Bu = temporary_fix_subsref(sys.b,1:sys.nx,sys.nw+1:sys.nw+sys.nu);
                sys.Dw = temporary_fix_subsref(sys.d,1:sys.nz,1:sys.nw); sys.Du = temporary_fix_subsref(sys.d,1:sys.nz,sys.nw+1:sys.nw+sys.nu);
            elseif sys.nu ~= 0 && sys.ny ~= 0 % output feedback case
                sys.Bw  = temporary_fix_subsref(sys.b,1:sys.nx,1:sys.nw);
                sys.Bu  = temporary_fix_subsref(sys.b,1:sys.nx,sys.nw+1:sys.nw+sys.nu);
                sys.Cz  = temporary_fix_subsref(sys.c,1:sys.nz,1:sys.nx);
                sys.Dzw = temporary_fix_subsref(sys.d,1:sys.nz,1:sys.nw);
                sys.Dzu = temporary_fix_subsref(sys.d,1:sys.nz,sys.nw+1:sys.nw+sys.nu);
                sys.Cy  = temporary_fix_subsref(sys.c,sys.nz+1:sys.nz+sys.ny,1:sys.nx);
                sys.Dyw = temporary_fix_subsref(sys.d,sys.nz+1:sys.nz+sys.ny,1:sys.nw);
                sys.Dyu = temporary_fix_subsref(sys.d,sys.nz+1:sys.nz+sys.ny,sys.nw+1:sys.nw+sys.nu);
            end
            
            function part = temporary_fix_subsref(M,rows,cols)
                % construct first row
                part = M(rows(1),cols(1)); % upper-left element
                for i = 2:length(cols) % concatenate columns to construct 1st row
                    part = [part, M(rows(1),cols(i))];
                end
                
                % construct row j and concatenate to previously constructed rows
                for j = 2:length(rows)
                    rowj = M(rows(j),cols(1)); % upper-left element
                    for i = 2:length(cols) % concatenate columns to construct 1st row
                        rowj = [rowj, M(rows(j),cols(i))];
                    end
                    part = [part; rowj];
                end
            end
        end
