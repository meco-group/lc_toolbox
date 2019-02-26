%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MULTI-OBJECTIVE OUTPUT FEEDBACK SYNTHESIS FOR LTI SYSTEMS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FOR THE MIMO LTI SYSTEMS
% dx = A x + Bw w + Bu u
% z = Cz x + Dzw w + Dzu u
% y = Cy x + Dyw w
% THE AIM IS TO DESING AN OUTPUT-FEEDBACK CONTROLLER OF THE FORM
% dxc = Ax xc + Bc y
% u = Cc xc + Dc y
% Hence, the resulting controller is of the same order as the generalized
% plant. Types of the performance objectives are H2/Hinf (also mixed)
% 
% INPUT ARGUMENTS
% sys:  instance of a state-space form class
% ny:   dimension of the measured output
% nu:   dimension of the control input
% opts: structure containing all the options required/extra options.
% polereg : specify pole region constraints.
%
% OUTPUT ARGUMENTS
% solution:     a structure providing the required solution/matrices and
% other important stuffs
% 
% Taranjitsingh Singh, Div. PMA, Dept. Mechanical Engineering, KU Leuven
% 2018
%
% References: 
%
% Carsten Schrer, Pascal Gahinet, and, Mahmoud Chilali
% Multiobjective Output-Feedback Control via LMI Optimization
% IEEE Transactions on Automatic Control, VOL. 42, No. 7, July 1997
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [K,sol_info] = LTI_unstructured_mix(sys,ny,nu,alpha,channel,polereg,varargin)

default = struct('spec', inf, 'tolerance',1e-6,'verbose',2, ...
    'scaling_obj',1,'solver','mosek','Dc',0);

if nargin == 7
    opts = mergestruct(varargin{1},default);
else
    opts = default;
end

t_start = cputime;

options = sdpsettings('solver',opts.solver,'verbose',opts.verbose,'mosek.MSK_DPAR_INTPNT_CO_TOL_REL_GAP',1e-6);
gen_sys = extract_generalized_plant(sys,nu,ny);

n_pspecs = length(channel);
A = gen_sys.A; nx = gen_sys.nx; Bu = gen_sys.Bu; Cy = gen_sys.Cy; 
[Bw,Cz,Dzw,Dzu,Dyw] = deal(cell(1,n_pspecs));
[nw,nz] = deal(zeros(1,n_pspecs));
        for j = 1:n_pspecs
            Bw{j}  = gen_sys.Bw(:,channel(j).In);
            Cz{j}  = gen_sys.Cz(channel(j).Out,:);   
            Dzw{j} = gen_sys.Dzw(channel(j).Out,channel(j).In); 
            Dzu{j} = gen_sys.Dzu(channel(j).Out,:); 
            Dyw{j} = gen_sys.Dyw(:,channel(j).In); 
            nw(j)  = size(Bw{j},2);
            nz(j)  = size(Cz{j},1);
            if sys.Ts == 0 && opts.spec == 2
                if any(Dzw(:)~=0)
                    disp('Continuous time, D-matrix closed-loop system nonzero => H2 performance infinite');
                    sol_info.feas      = 0;
                    sol_info.objective = [];
                    K = [];
                    return;
                end
            end
        end
        X = sdpvar(nx);
        Y = sdpvar(nx);
        Ac_hat = sdpvar(nx,nx,'full');
        Bc_hat = sdpvar(nx,ny,'full');
        Cc_hat = sdpvar(nu,nx,'full');
        if opts.Dc
        Dc_hat = sdpvar(nu,ny,'full');
        else
        Dc_hat = zeros(nu,ny);
        end
            for j = 1:n_pspecs
            gam2{j} = sdpvar(1);
            if channel(j).performance == 2
                W{j} = sdpvar(nz(j),nz(j),'full');
            end
        end
                Q = [X, eye(nx); eye(nx), Y];
        for j = 1:n_pspecs
            if channel(j).performance == inf                
                T11 = A*X + X*A' + Bu*Cc_hat + Cc_hat'*Bu';
                T21 = Ac_hat + (A + Bu*Dc_hat*Cy)';
                T22 = Y*A + A'*Y + Bc_hat*Cy + Cy'*Bc_hat';
                T31 = (Bw{j} + Bu*Dc_hat*Dyw{j})';
                T32 = (Y*Bw{j} + Bc_hat*Dyw{j})';
                T33 = -eye(nw(j));
                T41 = Cz{j}*X + Dzu{j}*Cc_hat;
                T42 = Cz{j} + Dzu{j}*Dc_hat*Cy;
                T43 = Dzw{j} + Dzu{j}*Dc_hat*Dyw{j};
                T44 = -gam2{j}*eye(nz(j));
                
                Term{j} = [T11, T21', T31', T41' ;
                           T21, T22 , T32', T42' ;
                           T31, T32 , T33 , T43' ; 
                           T41, T42 , T43 , T44 ];
            elseif channel(j).performance == 2
                T11 = A*X + X*A' + Bu*Cc_hat + Cc_hat'*Bu';
                T21 = Ac_hat + (A + Bu*Dc_hat*Cy)';
                T22 = Y*A + A'*Y + Bc_hat*Cy + Cy'*Bc_hat';
                T31 = (Bw{j} + Bu*Dc_hat*Dyw{j})';
                T32 = (Y*Bw{j} + Bc_hat*Dyw{j})';
                T33 = -eye(nw(j));
                
                Term1 = [T11, T21', T31' ;
                         T21, T22 , T32' ;
                         T31, T32 , T33 ];   
                        
                T11 = Q;     
                T21 = [Cz{j}*X + Dzu{j}*Cc_hat, Cz{j} + Dzu{j}*Dc_hat*Cy];
                T22 = W{j};   
                
                Term2 = [T11, T21' ; 
                         T21, T22 ];  
                     
                Term3 = trace(W{j}) - gam2{j};
                Term{j} = blkdiag(Term1,-Term2,Term3);
            end        
        end
Phi = [A*X+Bu*Cc_hat A;Ac_hat Y*A+Bc_hat*Cy];
comp_time_parsing = cputime - t_start;

% constraints
constraints = [Q >= 0];
for j = 1:n_pspecs
    constraints = [constraints, Term{j} <= 0];
    if alpha(j)==0
       constraints = [constraints, gam2{j} <= 1]; 
    end
end

n_pole_cons = length(polereg);
for j = 1:n_pole_cons; 
Pole{j} = kron(polereg{j}.L,Q) + kron(polereg{j}.M,Phi) + kron(polereg{j}.M',Phi');
constraints = [constraints, Pole{j} <= 1];
end

%objective 
objective = alpha(:)'*[gam2{:}]'*opts.scaling_obj;
%SolvetheSDP
sol = optimize(constraints,objective,options);
comp_time_SDP = cputime - t_start - comp_time_parsing;

%% extract solution info 

for j = 1:length(gam2)
    sol_info.gamma(j) = sqrt(value(gam2{j}));
end
sol_info.objective = value(objective);

X_ = value(X);
Y_ = value(Y);

A_ = value(Ac_hat);
B_ = value(Bc_hat);
C_ = value(Cc_hat);
D_ = value(Dc_hat);

[U,sig,V] = svd(eye(nx) - X_*Y_);

R = U*sqrt(sig);
S = V*sqrt(sig);

Dk = value(Dc_hat);
Ck = (value(Cc_hat) - Dk*Cy*X) * inv(R');
Bk = inv(S)*(value(Bc_hat) - Y*Bu*Dk);
Ak = inv(S)*(value(Ac_hat) - S*Bk*Cy*value(X) - value(Y)*Bu*Ck*R' - value(Y)*(A+Bu*Dk*Cy)*value(X))*inv(R');
K = SSmod(Ak,Bk,Ck,Dk);
K.Ts = sys.Ts;
sol_info.X = X_;
sol_info.Y = Y_;

if (all(eig(value(Q)) < 0))
error('conditions violated')
end

sol_info.comp_time_parsing = comp_time_parsing;
sol_info.comp_time_SDP     = comp_time_SDP;
sol_info.comp_time_total   = cputime - t_start;
             
end