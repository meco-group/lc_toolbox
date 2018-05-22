%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UNSTRUCTURED OUTPUT FEEDBACK SYNTHESIS FOR LPV SYSTEMS USING BSPLINES 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function should be able to handle multivariate case
% For the MIMO LPV system
%   dx =  A(a)x +  Bw(a)w +  Bu(a)u
%   z  = Cz(a)x + Dzw(a)w + Dzu(a)u
%   y  = Cy(a)x + Dyw(a)w 
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
% opts       : structure providing options
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
function [K,sol_info] = LPV_unstructured_OF(sys,ny,nu,param,varargin)

default = struct('spec',inf,'objective','L1','var_deg',1,...
    'var_knots',0,'relax_deg',0,'relax_mp',0,'tolerance',1e-6,...
    'verbose',0,'scaling_obj',1,'solver','mosek',...
    'controller_dependency','a');

if nargin == 5
    opts = mergestruct(varargin{1},default);
else
    opts = default;
end
switch sys.Ts
    case {0,{}}
        switch opts.controller_dependency
            case 'ada' % controller depending on parameter and its derivative
                [K,sol_info] = LPV_unstructured_OF_primal(sys,ny,nu,param,opts,0);
            case 'a'   % controller depending on parameter only
                [K1,Primal1] = LPV_unstructured_OF_primal(sys,ny,nu,param,opts,1);
                [K2,Primal2] = LPV_unstructured_OF_primal(sys,ny,nu,param,opts,2);
                if ((Primal1.objective <= Primal2.objective) && ~isinf(Primal1.objective))
                    sol_info = Primal1; K = K1;
                elseif Primal1.objective > Primal2.objective
                    sol_info = Primal2; K = K2;
                else
                    error('Sorry, SDP : INFEASIBLE problem, cannot go further')
                end
            otherwise
                error('Wrong selection of controller dependency, must be ''a'' or ''ada''.');
        end
    otherwise
        switch opts.controller_dependency
            case 'ada' % controller depending on parameter and its derivative
                [K,sol_info] = LPV_unstructured_OF_primal_DT(sys,ny,nu,param,opts,0);
            case 'a'   % controller depending on parameter only
                [K,sol_info] = LPV_unstructured_OF_primal_DT(sys,ny,nu,param,opts,1);
%                 [K2,Primal2] = LPV_unstructured_OF_primal_DT(sys,ny,nu,param,opts,2);
%                 if ((Primal1.objective <= Primal2.objective) && ~isinf(Primal1.objective))
%                     sol_info = Primal1; K = K1;
%                 elseif Primal1.objective > Primal2.objective
%                     sol_info = Primal2; K = K2;
%                 else
%                     error('Sorry, SDP : INFEASIBLE problem, cannot go further')
%                 end
%             otherwise
%                 error('Wrong selection of controller dependency, must be ''a'' or ''ada''.');
        end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [K,sol_info] = LPV_unstructured_OF_primal(sys,ny,nu,param,opts,dep)
    
check_optispline;
import splines.*;

t_start = cputime; 

options = sdpsettings('solver',opts.solver,'verbose',opts.verbose);
opti = OptiSplineYalmip();

% generalized plant matrices
gen_sys = extract_generalized_plant(sys,nu,ny);
A  = gen_sys.A;  Bw  = gen_sys.Bw;  Bu  = gen_sys.Bu;
Cz = gen_sys.Cz; Dzw = gen_sys.Dzw; Dzu = gen_sys.Dzu;
Cy = gen_sys.Cy; Dyw = gen_sys.Dyw; 
nx = gen_sys.nx; nz = gen_sys.nz; nw = gen_sys.nw;

% check to assure finite H2 norm
if sys.Ts == 0 && opts.spec == 2
    if isa(Dzw,'splines.Function')
        data = Dzw.coeff.data;
        if any(data(:)~=0)
            disp('Continuous time, D-matrix closed-loop system nonzero => H2 performance infinite');
            sol_info.feas      = 0;
            sol_info.objective = [];
            K = [];
            return;
        end
    elseif any(Dzw(:)~=0)
        disp('Continuous time, D-matrix closed-loop system nonzero => H2 performance infinite');
        sol_info.feas      = 0;
        sol_info.objective = []; 
        K = [];
        return;  
    end
end
% Theoretically speaking, the closed-loop D-matrix is given by
% Dcl = Dzw + Dzu*Dc*Dyw
% so in some cases Dcl can be zero if Dzw is nonzero. However, the latter
% case is omitted now for the sake of simplicity.
% In case of a strictly proper controller (i.e., Dc = 0), Dcl = 0 if and
% only if Dzw is zero.

%% Optimization variables
N = length(param);

% determine degree and knots for each basis
if     length(opts.var_deg) == 1; deg = opts.var_deg*ones(1,N);
elseif length(opts.var_deg) == N; deg = opts.var_deg;
else   error('opts.var_deg has incorrect dimension'); end
if     length(opts.var_knots) == 1; nknots = opts.var_knots*ones(1,N);
elseif length(opts.var_knots) == N; nknots = opts.var_knots;
else   error('opts.var_knots has incorrect dimension'); end

% construct TensorBasis
        for i = 1:N
            args{i} = param{i}.tensor_basis.arguments;
            B{i} = BSplineBasis([param{i}.basis.domain.min*ones(deg(i),1);linspace(param{i}.basis.domain.min,param{i}.basis.domain.max,nknots(i)+2)';param{i}.basis.domain.max*ones(deg(i),1)],deg(i));
        end
            TB = TensorBasis(B,args);
        
% construct LMI variables  
switch dep
    case 0 % X and Y dependent on constant and BRV parameters
        [X,dXdt] = build_Lyapmat(TB,param,nx,opti,'brv','c');
        [Y,dYdt] = build_Lyapmat(TB,param,nx,opti,'brv','c');
    case 1 % X depends on constant parameters only
        [X,dXdt] = build_Lyapmat(TB,param,nx,opti,'constant','c');
        [Y,dYdt] = build_Lyapmat(TB,param,nx,opti,'brv','c');
    case 2 % Y depends on constant parameters only
        [X,dXdt] = build_Lyapmat(TB,param,nx,opti,'brv','c');
        [Y,dYdt] = build_Lyapmat(TB,param,nx,opti,'constant','c');
    otherwise
        error ('something wrong in script or addition feature not implemented');
end
Ac_hat = opti.Function(TB,[nx,nx],'full');
Bc_hat = opti.Function(TB,[nx,ny],'full');
Cc_hat = opti.Function(TB,[nu,nx],'full');
Dc_hat = zeros(nu,ny);
switch opts.spec
    case inf
        switch opts.objective
            case 'wc'; gam2 = opti.variable(1);
            case 'L1'; gam2 = opti.Function(TB,[1,1],'symmetric');
            otherwise; error('PRIMAL SDP: invalid performance objective');
        end
    case 2
        Ny = opti.Function(TB,[nz,nz],'symmetric');
        switch opts.objective
            case 'wc'; gam2 = opti.variable(1);
            case 'L1'; gam2 = opti.Function(TB,[1,1],'symmetric');
            otherwise; error('PRIMAL SDP: invalid performance objective');
        end
end

%% Parameter-dependent LMIs
Q = [X, eye(gen_sys.nx); eye(gen_sys.nx), Y];

        switch opts.spec
            case 0 % stabilizing controller design
                        T11 = A*X + X*A' - dXdt + Bu*Cc_hat + Cc_hat'*Bu';
                        T21 = Ac_hat + (A + Bu*Dc_hat*Cy)';
                        T22 = Y*A + A'*Y + dYdt + Bc_hat*Cy + Cy'*Bc_hat';

                        Term = [T11, T21';
                                T21, T22];
            case inf % H-infinity controller design
                        T11 = A*X + X*A' - dXdt + Bu*Cc_hat + Cc_hat'*Bu';
                        T21 = Ac_hat + (A + Bu*Dc_hat*Cy)';
                        T22 = Y*A + A'*Y + dYdt + Bc_hat*Cy + Cy'*Bc_hat';
                        T31 = (Bw + Bu*Dc_hat*Dyw)';
                        T32 = (Y*Bw + Bc_hat*Dyw)';
                        T33 = -eye(nw);
                        T41 = Cz*X + Dzu*Cc_hat;
                        T42 = Cz + Dzu*Dc_hat*Cy;
                        T43 = Dzw + Dzu*Dc_hat*Dyw;
                        T44 = -gam2*eye(nz);

                        Term = [T11, T21', T31', T41' ;
                        T21, T22 , T32', T42' ;
                        T31, T32 , T33 , T43' ;
                        T41, T42 , T43 , T44 ];
            case 2 % H2 controller design

                        T11 = A*X + X*A' - dXdt + Bu*Cc_hat + Cc_hat'*Bu';
                        T21 = Ac_hat + (A + Bu*Dc_hat*Cy)';
                        T22 = Y*A + A'*Y + dYdt + Bc_hat*Cy + Cy'*Bc_hat';
                        T31 = (Bw + Bu*Dc_hat*Dyw)';
                        T32 = (Y*Bw + Bc_hat*Dyw)';
                        T33 = -eye(nw);  

                        Term1 = [T11, T21', T31' ;
                                 T21, T22 , T32' ;
                                 T31, T32 , T33 ];   

                        T11 = Q;     
                        T21 = [Cz*X + Dzu*Cc_hat, Cz + Dzu*Dc_hat*Cy];
                        T22 = Ny;           

                        Term2 = [T11, T21' ; 
                                 T21, T22 ];  

                        Term3 = trace(Ny) - gam2;

                        Term = blkdiag(Term1,-Term2,Term3);

                        % set Dc_hat = 0 to avoid constraint Dzw + Dzu*Dc_hat*Dyw == 0
                        % incorporate this later...
                        Dc_hat = zeros(nu,ny);
            Dc_hat = zeros(nu,ny);
        end     
%% LMI relaxations

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
    
%% solve SDP
constraints = {Q >= 0,Term <= 0};
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

% REMARK:
% the opti.solver command consumes a lot of time, while sol.solve() is fast

%% extract solution info
sol_info.coeffs_Q    = Q.coeff.dimension;
sol_info.coeffs_Term = Term.coeff.dimension;
[pr,~]               = checkset(opti.yalmip_constraints); % primal residual
if pr > -opts.tolerance
    sol_info.feasible  = 1;
    sol_info.objective = sol.value(objective);
    switch opts.spec
        case {2,inf}
            sol_info.gam2 = sol.value(gam2);
    end
    sol_info.X      = sol.value(X);
    sol_info.Y      = sol.value(Y);
    sol_info.dXdt   = sol.value(dXdt);
    sol_info.dYdt   = sol.value(dYdt);
    sol_info.Q      = sol.value(Q);
    sol_info.Ac_hat = sol.value(Ac_hat);
    sol_info.Bc_hat = sol.value(Bc_hat);
    sol_info.Cc_hat = sol.value(Cc_hat);
    sol_info.Dc_hat = sol.value(Dc_hat);
    sol_info.Term   = sol.value(Term);
else
    sol_info.feas      = 0;
    sol_info.objective = inf; 
    K = [];
    disp('PRIMAL SDP: infeasible problem'); return;
end

%% check solution
for j = 1:prod(sol_info.Q.coeff.dimension)
    min_eig(j) = min(eig(sol_info.Q.coeff.data(j))); 
end
sol_info.check_pos_Q = min(min_eig);
for j = 1:prod(sol_info.Term.coeff.dimension)
    max_eig(j) = max(eig(sol_info.Term.coeff.data(j))); 
end
sol_info.check_neg_Term = max(max_eig);    

%% CONTROLLER RECONTRUCTION
% The controller ss-matrices are reconstructed as follows:
%
% Dc = Dc_hat;
% Cc = (Cc_hat - Dc*Cy*X)*(M')^(-1);
% Bc = N^(-1)*(Bc_hat - Y*Bu*Dc);
% Ac = N^(-1)*(Ac_hat + Y*dXdt + N*dMdt' - Y*(A+Bu*Dc*Cy)*X - N*Bc*Cy*X - Y*Bu*Cc*M')*(M')^(-1);
%
% where N and M are solutions of the equation: I - X*Y = M*N'
%
% Note that the resulting controller depends on the parameter derivative!
%
% To obtain a 'practical' controller, run the above script with the
% one of the following changes to the optimization variables (such that 
% the term 
%     dXdt*Y + dMdt*N' = -(X*dYdt + M*dNdt') 
% in the reconstruction of Ac vanishes):
% 1) Take X independent of bRV parameters (hence dXdt = 0). Subsequently, select 
%        M = eye(nx) 
%        N = eye(nx) - Y*X 
%    for controller reconstruction. 
% 2) Take Y independent of bRV parameters (hence dYdt = 0). Subsequently, select
%        M = eye(nx) - X*Y
%        N = eye(nx)
%    for controller reconstruction.

% Extracting Lyapunov and slack variables
X_ = sol_info.X;
Y_ = sol_info.Y;
dX_dt = sol_info.dXdt;
dY_dt = sol_info.dYdt;
A_ = sol_info.Ac_hat;
B_ = sol_info.Bc_hat;
C_ = sol_info.Cc_hat;
D_ = sol_info.Dc_hat;

u = size(sys.B,2); u = u + ((-nu+1):0);
y = size(sys.C,1); y = y + ((-ny+1):0);

PCy_ = sys.C(y,:);
PBu_ = sys.B(:,u);
PA_  = sys.A;

if dep == 0     % controller depending on parameter and its derivative
    M_var = -X_*dY_dt;
    M11 = (X_*Y_)';
    M44 = zeros(nx);
elseif dep == 1 % controller depending on parameter only: X constant
    M_var = dX_dt*Y_;
    M11 = zeros(nx);
    M44 = Y_*X_;
elseif dep == 2 % controller depending on parameter only: Y constant
    M_var = -X_*dY_dt;
    M11 = (X_*Y_)';
    M44 = zeros(nx);
end

nx = size(A_,1);
Dt = D_;
Ct = (C_ - D_*PCy_*X_);
Bt = (B_ - Y_*PBu_*D_);
At = (A_ + M_var - Y_*(PA_+PBu_*D_*PCy_)*X_ - Bt*PCy_*X_ - Y_*PBu_*Ct);

nut = size(Bt,2);
nyt = size(Ct,1);

% Constructing LFT Matrices M,Nu and Nl                
M = {M11      , eye(nx)     , zeros(nx,nut), zeros(nx)    ;...
     zeros(nx), zeros(nx)   , zeros(nx,nut), eye(nx)      ;...
     Ct       , zeros(nyt,nx), Dt           , zeros(nyt,nx);...
     At       , zeros(nx)   , Bt           , M44         };
Nu = eye(nx); 
Nl = eye(nx);

param_new = [];
for i = 1:length(param)
    if (param{i}.bounded && dep == 0)
        param_new = [param_new,param{i}.add_parameter_derivative('c')];
    else
        param_new = [param_new,{param{i}}];
    end
end

% obtaining controller as an LFTmod object 
K = LFTmod(M,Nu,Nl,eye(nx),param_new);
K.Ts = sys.Ts;

% store computation times
sol_info.comp_time_parsing = comp_time_parsing;
sol_info.comp_time_SDP     = comp_time_SDP;
sol_info.comp_time_total   = cputime - t_start;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [K,sol_info] = LPV_unstructured_OF_primal_DT(sys,ny,nu,param,opts,dep)

% References :

% F. Amato, M. Mattei and A. Pironti.
% Gain scheduled control for discrete-time systems depending on bounded
% rate parameters.
% International Journal of Robust and Nonlinear Control (2005); 15:473-494

check_optispline;
import splines.*;

t_start = cputime; 

options = sdpsettings('solver',opts.solver,'verbose',opts.verbose);
opti = OptiSplineYalmip();

% generalized plant matrices
gen_sys = extract_generalized_plant(sys,nu,ny);
A  = gen_sys.A;  Bw  = gen_sys.Bw;  Bu  = gen_sys.Bu;
Cz = gen_sys.Cz; Dzw = gen_sys.Dzw; Dzu = gen_sys.Dzu;
Cy = gen_sys.Cy; Dyw = gen_sys.Dyw; 
nx = gen_sys.nx; nz = gen_sys.nz; nw = gen_sys.nw;

% Optimization variables
N = length(param);

% Taking care of parameter rates, if they are mentioned while defining
% parameter, see References.
for i = 1:N
    if any(abs(param{i}.rate) >= 1)
        error('You mentioned wrong bounds for rates of the any of the parameter for discrete time, it should be between -1 and 1, or do not put any');
    end
end

% determine degree and knots for each basis
if     length(opts.var_deg) == 1; deg = opts.var_deg*ones(1,N);
elseif length(opts.var_deg) == N; deg = opts.var_deg;
else   error('opts.var_deg has incorrect dimension'); end
if     length(opts.var_knots) == 1; nknots = opts.var_knots*ones(1,N);
elseif length(opts.var_knots) == N; nknots = opts.var_knots;
else   error('opts.var_knots has incorrect dimension'); end

        for i = 1:N
            args{i} = param{i}.tensor_basis.arguments;
            B{i} = BSplineBasis([param{i}.basis.domain.min*ones(deg(i),1);linspace(param{i}.basis.domain.min,param{i}.basis.domain.max,nknots(i)+2)';param{i}.basis.domain.max*ones(deg(i),1)],deg(i));
        end
            TB = TensorBasis(B,args);

% Construct LMI Variables
        switch dep
            case 0
                [X,dX] = build_Lyapmat(TB,param,nx,opti,'brv','d',deg,nknots);
                [Y,dY] = build_Lyapmat(TB,param,nx,opti,'brv','d',deg,nknots);
            case 1
                [X,dX] = build_Lyapmat(TB,param,nx,opti,'constant','d');
                [Y,dY] = build_Lyapmat(TB,param,nx,opti,'brv','d',deg,nknots);
            case 2
                [X,dX] = build_Lyapmat(TB,param,nx,opti,'brv','d',deg,nknots);
                [Y,dY] = build_Lyapmat(TB,param,nx,opti,'constant','d');
            otherwise
                error ('something wrong in script or addition feature not implemented');
        end
switch dep
    case {0,1}
        Ac_hat = opti.Function(X.tensor_basis,[nx,nx],'full');
    case 2
        Ac_hat = opti.Function(Y.tensor_basis,[nx,nx],'full');
end
Bc_hat = opti.Function(TB,[nx,ny],'full');
Cc_hat = opti.Function(TB,[nu,nx],'full');
Dc_hat = zeros(nu,ny);
switch opts.spec
    case inf
        switch opts.objective
            case 'wc'; gam2 = opti.variable(1);
            case 'L1'
                switch dep
                    case {0,1}
                        gam2 = opti.Function(Y.tensor_basis,[1,1],'symmetric');
                    case {2}
                        gam2 = opti.Function(X.tensor_basis,[1,1],'symmetric');
                end
            otherwise; error('PRIMAL SDP: invalid performance objective');
        end
    case 2
        W = opti.Function(Y.tensor_basis,[nz,nz],'symmetric');
        switch opts.objective
            case 'wc'; gam2 = opti.variable(1);
            case 'L1'
                switch dep
                    case {0,1}
                        gam2 = opti.Function(Y.tensor_basis,[1,1],'symmetric');
                    case {2}
                        gam2 = opti.Function(X.tensor_basis,[1,1],'symmetric');
                end
            otherwise; error('PRIMAL SDP: invalid performance objective');
        end
end
%% Parameter-dependent LMIs
Q = [X, eye(gen_sys.nx); eye(gen_sys.nx), Y];
switch opts.spec
    case 0 % stabilizing controller design
        
    T11 = dX;
    T12 = eye(nx);
    T13 = A*X+Bu*Cc_hat;
    T14 = A+Bu*Dc_hat*Cy;
    T22 = dY;
    T23 = Ac_hat;
    T24 = dY*A+Bc_hat*Cy;
    T33 = X;
    T34 = eye(nx);
    T44 = Y;
        
    Term = [T11 , T12 , T13 , T14;
            T12', T22 , T23 , T24;
            T13', T23', T33 , T34;
            T14', T24', T34', T44];
    Term = -Term;    
    case inf % H-infinity controller design
        
    T11 = dX;
    T12 = eye(nx);
    T13 = A*dX+Bu*Cc_hat;
    T14 = A+Bu*Dc_hat*Cy;
    T15 = Bw+Bu*Dc_hat*Dyw;
    T16 = zeros(nx,nz);
    T22 = dY;
    T23 = Ac_hat;
    T24 = dY*A+Bc_hat*Cy;
    T25 = dY*Bw+Bc_hat*Dyw;
    T26 = zeros(nx,nz);
    T33 = X;
    T34 = eye(nx);
    T35 = zeros(nx,nw);
    T36 = (Cz*dX+Dzu*Cc_hat)';
    T44 = Y;
    T45 = zeros(nx,nw);
    T46 = (Cz+Dzu*Dc_hat*Cy)';
    T55 = gam2*eye(nw);
    T56 = (Dzw+Dzu*Dc_hat*Dyw)';
    T66 = gam2*eye(nz);
    
    Term = [T11 , T12 , T13 , T14 , T15 , T16;
            T12', T22 , T23 , T24 , T25 , T26;
            T13', T23', T33 , T34 , T35 , T36;
            T14', T24', T34', T44 , T45 , T46;
            T15', T25', T35', T45', T55 , T56;
            T16', T26', T36', T46', T56', T66];
    Term = - Term;

    case 2 % H2 controller design                    
    T11 = dX;
    T12 = eye(nx);
    T13 = A*dX+Bu*Cc_hat;
    T14 = A+Bu*Dc_hat*Cy;
    T15 = Bw+Bu*Dc_hat*Dyw;
    T22 = dY;
    T23 = Ac_hat;
    T24 = dY*A+Bc_hat*Cy;
    T25 = dY*Bw+Bc_hat*Dyw;
    T33 = X;
    T34 = eye(nx);
    T35 = zeros(nx,nw);
    T44 = Y;
    T45 = zeros(nx,nw);
    T55 = gam2*eye(nw);
    
    Term1 = [T11 , T12 , T13 , T14 , T15;
            T12', T22 , T23 , T24 , T25;
            T13', T23', T33 , T34 , T35;
            T14', T24', T34', T44 , T45;
            T15', T25', T35', T45', T55];
    T11 = W;
    T12 = Cz*dX+Dzu*Cc_hat;
    T13 = Cz+Dzu*Dc_hat*Cy;
    T14 = Dzw+Dzu*Dc_hat*Dyw;
    T22 = X;
    T23 = eye(nx);
    T24 = zeros(nx,nw);
    T33 = Y;
    T34 = zeros(nx,nw);
    T44 = eye(nw);
    
    Term2 = [T11 , T12 , T13 , T14;
            T12', T22 , T23 , T24;
            T13', T23', T33 , T34;
            T14', T24', T34', T44];
                     
                Term3 = trace(W) - gam2;
                
                Term = blkdiag(-Term1,-Term2,Term3);

                % set Dc_hat = 0 to avoid constraint Dzw + Dzu*Dc_hat*Dyw == 0
                % incorporate this later...
                Dc_hat = zeros(nu,ny);
end

%% LMI relaxations

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
%% solve the sdp
constraints = {Q >= 0,Term <= 0};
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
%% extract solution info
sol_info.coeffs_Q    = Q.coeff.dimension;
sol_info.coeffs_Term = Term.coeff.dimension;
[pr,~]               = checkset(opti.yalmip_constraints); % primal residual
if pr > -opts.tolerance
    sol_info.feasible  = 1;
    sol_info.objective = sol.value(objective);
    switch opts.spec
        case {2,inf}
            sol_info.gam2 = sol.value(gam2);
    end
    sol_info.X      = sol.value(X);
    sol_info.Y      = sol.value(Y);
    sol_info.dX   = sol.value(dX);
    sol_info.dY   = sol.value(dY);
    sol_info.Q      = sol.value(Q);
    sol_info.Ac_hat = sol.value(Ac_hat); 
    sol_info.Bc_hat = sol.value(Bc_hat);
    sol_info.Cc_hat = sol.value(Cc_hat);
    sol_info.Dc_hat = sol.value(Dc_hat);
    sol_info.Term   = sol.value(Term);
else
    sol_info.feas      = 0;
    sol_info.objective = inf; 
    K = [];
    error('PRIMAL SDP: infeasible problem'); return;
end

%% check solution
% for j = 1:prod(sol_info.Q.coeff.dimension)
%     min_eig(j) = min(eig(sol_info.Q.coeff.data(j))); 
% end
% sol_info.check_pos_Q = min(min_eig);
% for j = 1:prod(sol_info.Term.coeff.dimension)
%     max_eig(j) = max(eig(sol_info.Term.coeff.data(j))); 
% end
% sol_info.check_neg_Term = max(max_eig);    

%% Controller Reconstruction
% Extracting Lyapunov and slack variables
X_ = sol_info.X;
Y_ = sol_info.Y;
dX_ = sol_info.dX;
dY_ = sol_info.dY;
A_ = sol_info.Ac_hat;
B_ = sol_info.Bc_hat;
C_ = sol_info.Cc_hat;
D_ = sol_info.Dc_hat;

u = size(sys.B,2); u = u + ((-nu+1):0);
y = size(sys.C,1); y = y + ((-ny+1):0);

PCy_ = sys.C(y,:);
PBu_ = sys.B(:,u);
PA_  = sys.A;
    if dep == 0     % controller depending on parameter and its derivative
        M11 = zeros(nx);
        M44 = dY_*X_;
    elseif dep == 1 % controller depending on parameter only: X constant
        M11 = zeros(nx);
        M44 = dY_*X_;
    elseif dep == 2 % controller depending on parameter only: Y constant
        M11 = dX_*Y_;
        M44 = zeros(nx);
    end
nx = size(A_,1);
Dt = D_;
Ct = (C_ - D_*PCy_*X_);
Bt = (B_ - Y_*PBu_*D_);
At = (A_ - dY_*(PA_+PBu_*D_*PCy_)*dX_ - Bt*PCy_*dX_ - dY_*PBu_*Ct);

nut = size(Bt,2);
nyt = size(Ct,1);

% Constructing LFT Matrices M,Nu and Nl                
M = {M11      , eye(nx)      , zeros(nx,nut), zeros(nx)    ;...
     zeros(nx), zeros(nx)    , zeros(nx,nut), eye(nx)      ;...
     Ct       , zeros(nyt,nx), Dt          , zeros(nyt,nx);...
     At       , zeros(nx)    , Bt           , M44         };
Nu = eye(nx); 
Nl = eye(nx);

param_new = [];
for i = 1:length(param)
    if (param{i}.bounded)
        param_new = [param_new,param{i}.add_parameter_derivative('d')];
    else
        param_new = [param_new,{param{i}}];
    end
end

% obtaining controller as an LFTmod object 
K = LFTmod(M,Nu,Nl,eye(nx),param_new);
K.Ts = sys.Ts;

% store computation times
sol_info.comp_time_parsing = comp_time_parsing;
sol_info.comp_time_SDP     = comp_time_SDP;
sol_info.comp_time_total   = cputime - t_start;
end
end