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

% function [Info] = compute_RO_mix(A,Bw,Bu,Cz,Cy,Dzw,Dzu,Dyw,sysKp,q,a,b,ch,opts,param)
%
% Compute a reduced-order multi-objective H2/Hinf controller for the 
% continuous or discrete-time LTI system 
%   dx = A  x + Bw  w + Bu  u
%   z  = Cz x + Dzw w + Dzu u
%   y  = Cy x + Dyw w + Dyu u
% Multiple H2 and/or Hinf performance criteria can be imposed: 
% A weighted average of optimization channels can be specified as objective
% function, and upper bounds on the H2/Hinf norms of each channel can be
% declared. 
%
% input: A,...,Dyw     -> generalized plant matrices
%        sysKp         -> controller of order p
%        q             -> controller order (q < p)
%        a             -> weighting vector for optimization channels. a(j) is 
%                         the weight for the H2 or Hinf norm of channel j. If 
%                         channel j is not an optimization channel, set a(j)=0.
%        b             -> vector with prespecified upper bounds. b(j) is an 
%                         upper bound for the H2 or Hinf norm of channel j. If 
%                         channel j is an optimization channel, set b(j)=0.
%        ch            -> Structure array 
%        ch.H2         -> Array specifying the H2 channels
%                         for example: channels.H2 = [1 2] means that
%                         an H2 specification is imposed on channels 1,2
%        ch.Hinf       -> Array specifying the Hinf channels
%        ch.In         -> Cell array containing selection matrices 
%                         corresponding to performance input w_j
%        ch.Out        -> Cell array containing selection matrices 
%                         corresponding to performance output z_j
%        opts          -> structure with options
%        opts.Dc       -> strictly proper controller design if opts.Dc = 0
%        param         -> struct with parameters
%        param.A22     -> tuning parameter in the LMI (default: 
%                         -I for continuous-time and 0 for discrete-time)
%        param.scaling -> LMI scaling parameter
%
% output: Info.feas    -> feasible (1) or not (0)
%         Info.cputime -> cpu time to solve the LMIs (seconds)
%         Info.LMIrows -> number of LMI rows
%         Info.LMIvars -> number of scalar variables
%         Info.P       -> Lyapunov matrix solutions
%         Info.X       -> LMI solution variable 
%         Info.K       -> stabilizing controller
%         Info.CL      -> closed-loop system
%         Info.mu      -> H2 upper bounds
%         Info.H2      -> H2 norms (computed a posteriori with norm.m)
%         Info.gam     -> Hinf upper bounds
%         Info.Hinf    -> Hinf norms (computed a posteriori with norm.m)
%         Info.stab    -> continuous-time : max(real(eig(Acl)))
%                         discrete-time   : max(abs(eig(Acl)))  
%
% References: 
%  G. Hilhorst, G. Pipeleers and J. Swevers.
%  Reduced-Order Multi-Objective H-infinity Control of an Overhead Crane 
%  Test Setup. In Proceedings of the 2013 CDC, Firenze, Italy. Dec 2013 
% 
% Gijs Hilhorst, Dept. of Mechanical Engineering, KU Leuven, 2013/2014
%
function [Info] = compute_RO_mix(A,Bw,Bu,Cz,Cy,Dzw,Dzu,Dyw,Dyu,sysKp,q,a,b,ch,opts,param)

% tolerance primal residual
tol = 1e-5;

% Total number of channels
nocH2   = length(ch.H2);
nocHinf = length(ch.Hinf);
noc     = nocH2 + nocHinf;

% convert controller to cell if necessary
if ~iscell(sysKp) % i.e. only one controller specified
    temp  = sysKp;
    sysKp = cell(1,noc);
    for j = 1:noc
        sysKp{j} = temp;
    end
end
Ts = sysKp{1}.Ts; % sampling time

% dimensions
nx = size(A,1); nu = size(Bu,2); ny = size(Cy,1); p = size(sysKp{1}.a,1);

% Error messages
if nargin < 14
    error('Input arguments incomplete')
elseif nargin == 14
    opts.Dc = 1;
    if Ts == 0
        A22 = -eye(p-q);
    else
        A22 = zeros(p-q,p-q);
    end
    scaling = 1;
elseif nargin == 15
    if ~isfield(opts,'Dc')
        opts.Dc = 1;
    end
    if Ts == 0
        A22 = -eye(p-q);
    else
        A22 = zeros(p-q,p-q);
    end
    scaling = 1;
elseif nargin == 16
    if ~isfield(opts,'Dc')
        opts.Dc = 1;
    end
    if isfield(param,'A22')
        A22 = param.A22;
    else
        if Ts == 0
            A22 = -eye(p-q);
        else
            A22 = zeros(p-q,p-q);
        end
    end
    if isfield(param,'scaling')
        scaling = param.scaling;
    else
        scaling = 1;
    end
elseif nargin > 15
    error('Too many input arguments')
end

% generalized plant
sysP = ss(A,[Bw,Bu],[Cz;Cy],[Dzw,Dzu;Dyw,Dyu],Ts);

% System matrices and dimensions for specific channels
[BBw,CCz,DDzw,DDzu,DDyw,nnw,nnz,sysPP] = deal(cell(1,noc));
for j=1:noc
    BBw{j}   = Bw*ch.In{j};
    CCz{j}   = ch.Out{j}*Cz;
    DDzw{j}  = ch.Out{j}*Dzw*ch.In{j};
    DDzu{j}  = ch.Out{j}*Dzu;
    DDyw{j}  = Dyw*ch.In{j};
    nnw{j}   = size(BBw{j},2);
    nnz{j}   = size(CCz{j},1);
    sysPP{j} = ss(A,[BBw{j},Bu],[CCz{j};Cy],[DDzw{j},DDzu{j};DDyw{j},Dyu],Ts);
end

% new LMI system
LMIs = [];

% controller parameters order p
Kp = cell(1,noc);
for j = 1:noc
    Ac = sysKp{j}.a; Bc = sysKp{j}.b;
    Cc = sysKp{j}.c; Dc = sysKp{j}.d;
    Kp{j} = [Ac,Bc;Cc,Dc];
end

% closed-loop interconnections sysPP{j} and sysKp{j}
[Acl,Bcl,Ccl,Dcl] = deal(cell(1,noc));
for j = 1:noc
    sysCLp = lft(sysPP{j},sysKp{j},nu,ny);
    Acl{j} = sysCLp.a; Bcl{j} = sysCLp.b;
    Ccl{j} = sysCLp.c; Dcl{j} = sysCLp.d;
end

% augmented matrices
[Bwa,Cza,Dzwa,Dzua,Dywa] = deal(cell(1,noc));
for j = 1:noc
    [Aa,Bwa{j},Bua,Cza{j},Cya,Dzwa{j},Dzua{j},Dywa{j}] = aug_model(A,BBw{j},Bu,CCz{j},Cy,DDzw{j},DDzu{j},DDyw{j},p);
end

% Scaling matrices
U = [eye(q)     ,zeros(q,p-q) ,zeros(q,nu);
     zeros(nu,q),zeros(nu,p-q),eye(nu)   ];
V = [eye(q)     ,zeros(q,p-q) ,zeros(q,ny);
     zeros(ny,q),zeros(ny,p-q),eye(ny)   ];

% LMI variables
X11 = sdpvar(q ,q ,'full');
X13 = sdpvar(q ,nu,'full');
X31 = sdpvar(nu,q ,'full');
X33 = sdpvar(nu,nu,'full');
[P,X12,X22,X32,X] = deal(cell(1,noc));
for j = 1:noc
    P{j}   = sdpvar(nx+p,nx+p,'symmetric');
    X12{j} = sdpvar(q  ,p-q,'full');
    X22{j} = sdpvar(p-q,p-q,'full');
    X32{j} = sdpvar(nu ,p-q,'full');
    X{j}   = [X11         ,X12{j},X13          ;
              zeros(p-q,q),X22{j},zeros(p-q,nu);
              X31         ,X32{j},X33         ];
end
M   = [zeros(q,q)  ,zeros(q,p-q) ,zeros(q,ny)  ;
       zeros(p-q,q),A22          ,zeros(p-q,ny);
       zeros(nu,q) ,zeros(nu,p-q),zeros(nu,ny)];
gam = cell(1,length(ch.Hinf));
for j = ch.Hinf
    gam{j} = sdpvar(1); % Hinf bounds
end
[W,mu] = deal(cell(1,length(ch.H2)));
for j = ch.H2
    W{j}  = sdpvar(nnz{j},nnz{j},'symmetric');
    mu{j} = sdpvar(1);  % H2 bounds
end
Kq_hat = sdpvar(q+nu,p+ny,'full');
if opts.Dc == 0 % strictly proper controller
   Kq_hat(q+1:end,q+1:end) = 0;
   X31 = zeros(nu,q);
end

Info.LMIrows = 0;
if Ts == 0 % continuous-time
    
    % Hinf LMIs
    for j = ch.Hinf   
        T11 = Acl{j}'*P{j} + P{j}*Acl{j};
        T12 = P{j}*Bcl{j};
        T13 = Ccl{j}';
        T14 = P{j}*Bua + Cya'*(U'*Kq_hat+X{j}*M-X{j}*Kp{j})';
        T22 = -gam{j}*eye(nnw{j});
        T23 = Dcl{j}';
        T24 = Dywa{j}'*(U'*Kq_hat+X{j}*M-X{j}*Kp{j})';
        T33 = -gam{j}*eye(nnz{j});
        T34 = Dzua{j};
        T44 = -X{j}-X{j}';
        
        Term = [T11 , T12 , T13 , T14;
                T12', T22 , T23 , T24;
                T13', T23', T33 , T34;
                T14', T24', T34', T44];
        
        sc1 = eye(nx+p+nnw{j}+nnz{j}); sc2 = scaling*eye(p+nu);
        Term = blkdiag(sc1,sc2)'*Term*blkdiag(sc1,sc2);   
        LMIs = [LMIs ,Term < 0];
        Info.LMIrows = Info.LMIrows + nx + p + nnw{j} + nnz{j} + p + nu;
    end
    
    % H2 LMIs 
    for j = ch.H2   
        T11 = Acl{j}'*P{j} + P{j}*Acl{j};
        T12 = P{j}*Bcl{j};
        T13 = P{j}*Bua + Cya'*(U'*Kq_hat+X{j}*M-X{j}*Kp{j})';
        T22 = -mu{j}*eye(nnw{j});
        T23 = Dywa{j}'*(U'*Kq_hat+X{j}*M-X{j}*Kp{j})';
        T33 = -X{j}-X{j}';

        Term = [T11 , T12 , T13;
                T12', T22 , T23;
                T13', T23', T33];

        sc1 = eye(nx+p+nnw{j}); sc2 = scaling*eye(p+nu);
        Term = blkdiag(sc1,sc2)'*Term*blkdiag(sc1,sc2);   
        LMIs = [LMIs ,Term < 0];
        Info.LMIrows = Info.LMIrows + nx + p + nnw{j} + p + nu; 

        T11 = W{j};
        T12 = Ccl{j};
        T13 = Dzua{j};
        T22 = P{j};
        T23 = -Cya'*(U'*Kq_hat+X{j}*M-X{j}*Kp{j})';
        T33 = X{j}+X{j}';

        Term = [T11 , T12 , T13;
                T12', T22 , T23;
                T13', T23', T33];
        
        sc1 = eye(nnz{j}+nx+p); sc2 = scaling*eye(p+nu);
        Term = blkdiag(sc1,sc2)'*Term*blkdiag(sc1,sc2);   
        LMIs = [LMIs ,Term > 0];
        Info.LMIrows = Info.LMIrows + nnz{j} + nx + p + p + nu;

        LMIs = [LMIs ,trace(W{j}) < mu{j}];
    end
else % discrete-time
    
    % H-infinity LMIs 
    for j = ch.Hinf   
        T11 = P{j};
        T12 = P{j}*Acl{j};
        T13 = P{j}*Bcl{j}; 
        T14 = zeros(nx+p,nnz{j});
        T15 = P{j}*Bua;
        T22 = P{j};
        T23 = zeros(nx+p,nnw{j});
        T24 = Ccl{j}';
        T25 = Cya'*(U'*Kq_hat+X{j}*M-X{j}*Kp{j})';
        T33 = gam{j}*eye(nnw{j});
        T34 = Dcl{j}';
        T35 = Dywa{j}'*(U'*Kq_hat+X{j}*M-X{j}*Kp{j})';
        T44 = gam{j}*eye(nnz{j});
        T45 = Dzua{j};
        T55 = -X{j}-X{j}';

        Term = [T11 , T12 , T13 , T14 , T15;
                T12', T22 , T23 , T24 , T25;
                T13', T23', T33 , T34 , T35;
                T14', T24', T34', T44 , T45;
                T15', T25', T35', T45', T55];

        sc1 = eye(2*nx+2*p+nnw{j}+nnz{j}); sc2 = scaling*eye(p+nu);
        Term = blkdiag(sc1,sc2)'*Term*blkdiag(sc1,sc2);   
        LMIs = [LMIs, Term > 0];
        Info.LMIrows = Info.LMIrows + 2*(nx+p) + nnw{j} + nnz{j} + p + nu;
    end
    
    % H2 LMIs
    for j = ch.H2
        T11 = P{j};
        T12 = P{j}*Acl{j};
        T13 = P{j}*Bcl{j}; 
        T14 = P{j}*Bua;
        T22 = P{j};
        T23 = zeros(nx+p,nnw{j});
        T24 = Cya'*(U'*Kq_hat+X{j}*M-X{j}*Kp{j})';
        T33 = mu{j}*eye(nnw{j});
        T34 = Dywa{j}'*(U'*Kq_hat+X{j}*M-X{j}*Kp{j})';
        T44 = -X{j}-X{j}';

        Term = [T11 , T12 , T13 , T14;
                T12', T22 , T23 , T24;
                T13', T23', T33 , T34;
                T14', T24', T34', T44];

        sc1 = eye(2*nx+2*p+nnw{j}); sc2 = scaling*eye(p+nu);
        Term = blkdiag(sc1,sc2)'*Term*blkdiag(sc1,sc2);   
        LMIs = [LMIs, Term > 0];
        Info.LMIrows = Info.LMIrows + 2*(nx+p) + nnw{j} + p + nu;

        T11 = W{j};
        T12 = Ccl{j};
        T13 = Dcl{j}; 
        T14 = Dzua{j};
        T22 = P{j};
        T23 = zeros(nx+p,nnw{j});
        T24 = Cya'*(U'*Kq_hat+X{j}*M-X{j}*Kp{j})';
        T33 = mu{j}*eye(nnw{j});
        T34 = Dywa{j}'*(U'*Kq_hat+X{j}*M-X{j}*Kp{j})';
        T44 = -X{j}-X{j}';

        Term = [T11 , T12 , T13 , T14;
                T12', T22 , T23 , T24;
                T13', T23', T33 , T34;
                T14', T24', T34', T44];

        sc1 = eye(nnz{j}+nx+p+nnw{j}); sc2 = scaling*eye(p+nu);
        Term = blkdiag(sc1,sc2)'*Term*blkdiag(sc1,sc2);   
        LMIs = [LMIs, Term > 0];
        Info.LMIrows = Info.LMIrows + nnz{j} + nx + p + nnw{j} + p + nu;

        LMIs = [LMIs, trace(W{j}) < mu{j}];
    end
end

% constraints
for j = ch.Hinf 
    if b(j) > 0 % Hinf bound channel j
        LMIs = [LMIs, gam{j} < b(j)];
    end
end
for j = ch.H2 
    if b(j) > 0 % H2 bound channel j
        LMIs = [LMIs, mu{j} < b(j)];
    end
end
Info.LMIvars = size(getvariables(LMIs),2);

% objective function: a*[gam{j};mu{j}]
obj = 0;
for j = ch.Hinf
    obj = obj + a(j)*gam{j};
end
for j = ch.H2
    obj = obj + a(j)*mu{j};
end

% Solve the SDP
sol = solvesdp(LMIs,obj,sdpsettings('verbose',0,'solver','mosek'));
Info.cpusec = sol.solvertime;
p_res = checkset(LMIs);

% extract solution when it exists
Info.feas = 0;
if p_res > -tol
    Info.feas = 1;
  
    % solution variables
    for j = 1:noc
        Info.P{j} = double(P{j});
        Info.X{j} = double(X{j});
    end
    Kq_hat = double(Kq_hat);
    X11 = double(X11); X13 = double(X13);
    X31 = double(X31); X33 = double(X33);
      
    % back-transformation of controller variables
    Kq = [X11, X13; X31, X33]\Kq_hat*V';
    Ac = Kq(1:q,1:q);     Bc = Kq(1:q,q+1:end);
    Cc = Kq(q+1:end,1:q); Dc = Kq(q+1:end,q+1:end);

    Info.Kq = ss(Ac,Bc,Cc,Dc,Ts);      % controller  
    Info.CL = lft(sysP,Info.Kq,nu,ny); % closed-loop system
    for j = ch.Hinf % Hinf upper bounds
        Info.gam(j) = double(gam{j});
    end
    for j = ch.H2   % H2 upper bounds
        Info.mu(j) = double(mu{j});
    end  
    sys.CLch = cell(1,noc);
    for j = ch.Hinf % a posteriori computed Hinf norms
        A = Info.CL.a;           B = Info.CL.b*ch.In{j};
        C = ch.Out{j}*Info.CL.c; D = ch.Out{j}*Info.CL.d*ch.In{j};
        sys.CLch{j} = ss(A,B,C,D,Ts);
        Info.Hinf(j) = norm(sys.CLch{j},inf);
    end
    for j = ch.H2   % a posteriori computed H2 norms
        A = Info.CL.a;           B = Info.CL.b*ch.In{j};
        C = ch.Out{j}*Info.CL.c; D = ch.Out{j}*Info.CL.d*ch.In{j};
        sys.CLch{j} = ss(A,B,C,[],Ts);
        Info.H2(j) = norm(sys.CLch{j});        
    end
    if Ts == 0
        Info.stab = max(real(eig(Info.CL.a)));
    else
        Info.stab = max(abs(eig(Info.CL.a)));
    end
end
 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Aa,Bwa,Bua,Cza,Cya,Dzwa,Dzua,Dywa] = aug_model(A,Bw,Bu,Cz,Cy,Dzw,Dzu,Dyw,nc)

% Dimensions
nx = size(A,1);
nw = size(Bw,2);
nu = size(Bu,2);
nz = size(Cz,1);
ny = size(Cy,1);

% Augmented system
Aa   = [A,zeros(nx,nc);zeros(nc,nx),zeros(nc,nc)];
Bwa  = [Bw;zeros(nc,nw)];
Bua  = [zeros(nx,nc),Bu;eye(nc),zeros(nc,nu)];
Cza  = [Cz,zeros(nz,nc)];
Cya  = [zeros(nc,nx),eye(nc);Cy,zeros(ny,nc)];
Dzwa = Dzw;
Dzua = [zeros(nz,nc),Dzu];
Dywa = [zeros(nc,nw);Dyw]; 

end

