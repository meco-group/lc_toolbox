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

% function [Info] = compute_FO_mix(A,Bw,Bu,Cz,Cy,Dzw,Dzu,Dyw,Dyu,a,b,ch,opts)
%
% Compute a full-order multi-objective H2/Hinf controller for the 
% continuous or discrete-time LTI system 
%   dx = A  x + Bw  w + Bu  u
%   z  = Cz x + Dzw w + Dzu u
%   y  = Cy x + Dyw w + Dyu u
% Multiple H2 and/or Hinf performance criteria can be imposed: 
% A weighted average of optimization channels can be specified as objective
% function, and upper bounds on the H2/Hinf norms of each channel can be
% declared. Both the Lyapunov shaping paradigm and the less conservative 
% G-shaping paradigm are implemented.
%
% input: A,...,Dyu   -> generalized plant matrices
%        a           -> weighting vector for optimization channels. a(j) is 
%                       the weight for the H2 or Hinf norm of channel j. If 
%                       channel j is not an optimization channel, set a(j)=0.
%        b           -> vector with prespecified upper bounds. b(j) is an 
%                       upper bound for the H2 or Hinf norm of channel j. 
%                       If channel j is an optimization channel, set b(j)=0.
%        ch          -> Structure array specifying performance channels 
%        ch.H2       -> Array specifying the H2 channels
%                       for example: channels.H2 = [1 2] means that
%                       an H2 specification is imposed on channels 1,2
%        ch.Hinf     -> Array specifying the Hinf channels
%        ch.In       -> Cell array containing selection matrices 
%                       corresponding to performance input w_j
%        ch.Out      -> Cell array containing selection matrices 
%                       corresponding to performance output z_j
%        opts        -> Struct with options
%        opts.Dc     -> strictly proper controller (0), else Dc~=0 
%        opts.method -> Lyapunov shaping 'P' or G-shaping 'G'
%
% output: Info.feas    -> feasible (1) or not (0)
%         Info.cputime -> cpu time to solve the LMIs (seconds)
%         Info.LMIrows -> number of LMI rows
%         Info.LMIvars -> number of scalar variables
%         Info.P       -> cell with Lyapunov matrix solutions
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
%
%   M.C. de Oliveira, J.C. Geromel and J. Bernussou.
%   Extended H2 and Hinf norm characterizations and controller
%   parametrizations for discrete-time systems (2002).
%   International Journal of Control. 2002;75(9):666-679
%
%   C.W. Scherer, P. Gahinet and M. Chilali.
%   Multiobjective Output-Feedback Control via LMI Optimization (1997).
%   IEEE Transactions on Automatic Control. 1997;42(7):896-911   
%
%   S. Boyd, L. El Ghaoui, E. Feron and V. Balakrishnan.
%   Linear Matrix Inequalities in System and Control Theory
%   vol. 15 of Studies in Applied Mathematics. Philadelphia: SIAM; 1994
%
% Gijs Hilhorst, Dept. of Mechanical Engineering, KU Leuven, 2013/2014
%
function [Info] = compute_FO_mix(A,Bw,Bu,Cz,Cy,Dzw,Dzu,Dyw,Dyu,a,b,ch,opts)

% Default options
if nargin < 12
    error('Input arguments incomplete')
elseif nargin == 12
    opts.Ts = 0;
    opts.Dc = 1;
    opts.method = 'P';
    opts.verbose = 1;
    opts.solver = 'mosek';
elseif nargin == 13
    if ~isfield(opts,'verbose'), opts.verbose = 1; end
    if ~isfield(opts,'solver'), opts.solver = 'mosek'; end
elseif nargin > 13
    error('Too many input arguments')
end

% Solve the controller design problem
if opts.method == 'P'     % Lyapunov shaping
    [Info] = Pshaping(A,Bw,Bu,Cz,Cy,Dzw,Dzu,Dyw,Dyu,a,b,ch,opts);
elseif opts.method == 'G' % G-shaping
    if opts.Ts == 0;
        error('G-shaping only for discrete-time models')
    end
    [Info] = Gshaping(A,Bw,Bu,Cz,Cy,Dzw,Dzu,Dyw,Dyu,a,b,ch,opts);
else
    error('Invalid solution method selected. Choose Lyapunov shaping or G-shaping');
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Info] = Pshaping(A,Bw,Bu,Cz,Cy,Dzw,Dzu,Dyw,Dyu,a,b,ch,opts)

tol     = 1e-5; % tolerance primal residual
tolsing = 1e-9; % tolerance algebraic loop

% dimensions
nx = size(A,1); nu = size(Bu,2); ny = size(Cy,1);

% generalized plant
sys.P = ss(A,[Bw,Bu],[Cz;Cy],[Dzw,Dzu;Dyw,Dyu],opts.Ts);

% Total number of channels
nocH2   = length(ch.H2);
nocHinf = length(ch.Hinf);
noc     = nocH2 + nocHinf;

% System matrices and dimensions for specific channels
[BBw,CCz,DDzw,DDzu,DDyw,nnw,nnz] = deal(cell(1,noc));
for j=1:noc
    BBw{j}  = Bw*ch.In{j};
    CCz{j}  = ch.Out{j}*Cz;
    DDzw{j} = ch.Out{j}*Dzw*ch.In{j};
    DDzu{j} = ch.Out{j}*Dzu;
    DDyw{j} = Dyw*ch.In{j};
    nnw{j}  = size(BBw{j},2);
    nnz{j}  = size(CCz{j},1);
end

% new LMI system
LMIs = [];
   
% LMI variables
X = sdpvar(nx,nx,'symmetric'); % Lyapunov matrix variables
Y = sdpvar(nx,nx,'symmetric');
Q = sdpvar(nx,nx,'full');      % transformed controller variables
F = sdpvar(nx,ny,'full');
L = sdpvar(nu,nx,'full');
R = sdpvar(nu,ny,'full');
if opts.Dc == 0;               % strictly proper controller
    R = zeros(nu,ny);
end
gam = cell(1,length(ch.Hinf));
for j = ch.Hinf
    gam{j} = sdpvar(1); % Hinf bounds
end
[W,mu] = deal(cell(1,length(ch.H2)));
for j = ch.H2
    W{j}  = sdpvar(nnz{j},nnz{j},'symmetric');
    mu{j} = sdpvar(1);  % H2 bounds
end

Info.LMIrows = 0;
if opts.Ts == 0 % continuous-time
    % Lyapunov matrix positive definite
    LMIs = [LMIs, [X,eye(nx);eye(nx),Y] > 0];
    Info.LMIrows = Info.LMIrows + 2*nx;
    
    % H-infinity LMIs 
    for j = ch.Hinf      
        T11 = A*X+X*A'+Bu*L+L'*Bu';
        T12 = Q'+A+Bu*R*Cy;
        T13 = BBw{j}+Bu*R*DDyw{j};
        T14 = (CCz{j}*X+DDzu{j}*L)';
        T22 = A'*Y+Y*A+F*Cy+Cy'*F';
        T23 = Y*BBw{j}+F*DDyw{j};
        T24 = (CCz{j}+DDzu{j}*R*Cy)';
        T33 = -gam{j}*eye(nnw{j});
        T34 = (DDzw{j}+DDzu{j}*R*DDyw{j})';
        T44 = -gam{j}*eye(nnz{j});
        
        Term = [T11 , T12 , T13 , T14;
                T12', T22 , T23 , T24;
                T13', T23', T33 , T34;
                T14', T24', T34', T44];
                   
        LMIs = [LMIs, Term <= 0];
        Info.LMIrows = Info.LMIrows + 2*nx + nnw{j} + nnz{j};
    end
    
    % H2 LMIs
    for j = ch.H2
        T11 = A*X+X*A'+Bu*L+L'*Bu';
        T12 = Q'+A+Bu*R*Cy;
        T13 = BBw{j}+Bu*R*DDyw{j};
        T22 = A'*Y+Y*A+F*Cy+Cy'*F';
        T23 = Y*BBw{j}+F*DDyw{j};
        T33 = -mu{j}*eye(nnw{j});

        Term = [T11 , T12 , T13;
                T12', T22 , T23;
                T13', T23', T33];  
        
        LMIs = [LMIs, Term <= 0];
        Info.LMIrows = Info.LMIrows + 2*nx + nnw{j}; 
           
        T11 = X;
        T12 = eye(nx);
        T13 = (CCz{j}*X+DDzu{j}*L)';
        T22 = Y;
        T23 = (CCz{j}+DDzu{j}*R*Cy)';
        T33 = W{j};

        Term = [T11 , T12 , T13;
                T12', T22 , T23;
                T13', T23', T33];  

        LMIs = [LMIs, Term >= 0];
        Info.LMIrows = Info.LMIrows + 2*nx + nnz{j};

        LMIs = LMIs + [trace(W{j}) <= mu{j}];       
        if opts.Dc ~= 0
            if norm(DDzu{j}) > 0 
                if norm(DDyw{j}) > 0 
                    LMIs = [LMIs, DDzw{j}+DDzu{j}*R*DDyw{j} == 0];
                end
            end
        end
    end
else % discrete-time
    
    % H-infinity LMIs 
    for j = ch.Hinf   
        T11 = X;
        T12 = eye(nx);
        T13 = A*X+Bu*L;
        T14 = A+Bu*R*Cy;
        T15 = BBw{j}+Bu*R*DDyw{j};
        T16 = zeros(nx,nnz{j});
        T22 = Y;
        T23 = Q;
        T24 = Y*A+F*Cy;
        T25 = Y*BBw{j}+F*DDyw{j};
        T26 = zeros(nx,nnz{j});
        T33 = X;
        T34 = eye(nx);
        T35 = zeros(nx,nnw{j});
        T36 = (CCz{j}*X+DDzu{j}*L)';
        T44 = Y;
        T45 = zeros(nx,nnw{j});
        T46 = (CCz{j}+DDzu{j}*R*Cy)';
        T55 = gam{j}*eye(nnw{j});
        T56 = (DDzw{j}+DDzu{j}*R*DDyw{j})';
        T66 = gam{j}*eye(nnz{j});

        Term = [T11 , T12 , T13 , T14 , T15 , T16;
                T12', T22 , T23 , T24 , T25 , T26;
                T13', T23', T33 , T34 , T35 , T36;
                T14', T24', T34', T44 , T45 , T46;
                T15', T25', T35', T45', T55 , T56;
                T16', T26', T36', T46', T56', T66];
                   
        LMIs = LMIs + [Term >= 0];
        Info.LMIrows = Info.LMIrows + 4*nx + nnw{j} + nnz{j};
    end
    
    % H2 LMIs
    for j = ch.H2
        T11 = X;
        T12 = eye(nx);
        T13 = A*X+Bu*L;
        T14 = A+Bu*R*Cy;
        T15 = BBw{j}+Bu*R*DDyw{j};
        T22 = Y;
        T23 = Q;
        T24 = Y*A+F*Cy;
        T25 = Y*BBw{j}+F*DDyw{j};
        T33 = X;
        T34 = eye(nx);
        T35 = zeros(nx,nnw{j});
        T44 = Y;
        T45 = zeros(nx,nnw{j});
        T55 = mu{j}*eye(nnw{j});

        Term = [T11 , T12 , T13 , T14 , T15;
                T12', T22 , T23 , T24 , T25;
                T13', T23', T33 , T34 , T35;
                T14', T24', T34', T44 , T45;
                T15', T25', T35', T45', T55];
    
        LMIs = LMIs + [Term >= 0];
        Info.LMIrows = Info.LMIrows + 4*nx + nnw{j};
           
        T11 = W{j};
        T12 = CCz{j}*X+DDzu{j}*L;
        T13 = CCz{j}+DDzu{j}*R*Cy;
        T14 = DDzw{j}+DDzu{j}*R*DDyw{j};
        T22 = X;
        T23 = eye(nx);
        T24 = zeros(nx,nnw{j});
        T33 = Y;
        T34 = zeros(nx,nnw{j});
        T44 = mu{j}*eye(nnw{j});

        Term = [T11 , T12 , T13 , T14;
                T12', T22 , T23 , T24;
                T13', T23', T33 , T34;
                T14', T24', T34', T44];
        
        LMIs = LMIs + [Term >= 0];
        Info.LMIrows = Info.LMIrows + nnz{j} + 2*nx + nnw{j};
            
        LMIs = LMIs + [trace(W{j}) <= mu{j}];
    end
end

% constraints
for j = ch.Hinf 
    if b(j) > 0 % Hinf bound channel j
        LMIs = LMIs + [gam{j} <= b(j)];
    end
end
for j = ch.H2 
    if b(j) > 0 % H2 bound channel j
        LMIs = LMIs + [mu{j} <= b(j)];
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
sol = solvesdp(LMIs,obj,sdpsettings('verbose',opts.verbose,'solver',opts.solver));
Info.cpusec = sol.solvertime;
p_res = checkset(LMIs);

% extract solution when it exists
Info.feas = 0;
if p_res > -tol
    Info.feas = 1;
    
    % solution variables
    X = double(X); 
    Y = double(Y);
    Q = double(Q); 
    F = double(F); 
    L = double(L); 
    if opts.Dc ~= 0
        R = double(R);
    end
    Info.P = [X,eye(nx);eye(nx),Y];
    
    % back-transformation of controller variables
    [Ac,Bc,Cc,Dc] = controller_transf(A,Bu,Cy,X,Y,Q,F,L,R,eye(nx)); 
  
    % compensation for Dyu ~= 0
    if norm(Dyu,1) > 0
        if norm(Dc,1) > 0
            Myuc = eye(ny)+Dyu*Dc; 
            Mcyu = eye(nu)+Dc*Dyu;
            SV = svd(Myuc);
            if min(SV) < sqrt(tolsing)
                error('Algebraic loop');
            else
                Cc = Mcyu\Cc;
                Dc = Mcyu\Dc;
                Ac = Ac-Bc*Dyu*Cc;
                Bc = Bc/Myuc;
            end
        else
            Ac = Ac-Bc*Dyu*Cc;
        end
    end
    Info.K = ss(Ac,Bc,Cc,Dc,opts.Ts);  % controller
    Info.CL = lft(sys.P,Info.K,nu,ny); % closed-loop system    
    for j = ch.Hinf % Hinf upper bounds
        Info.gam(j) = double(gam{j});
    end
    for j = ch.H2   % H2 upper bounds
        Info.mu(j) = double(mu{j});
    end  
    sys.CLch = cell(1,noc);
    for j = ch.Hinf % a posteriori computed Hinf norms
        A = Info.CL.a;
        B = Info.CL.b*ch.In{j};
        C = ch.Out{j}*Info.CL.c; 
        D = ch.Out{j}*Info.CL.d*ch.In{j};
        sys.CLch{j} = ss(A,B,C,D,opts.Ts);
        Info.Hinf(j) = norm(sys.CLch{j},inf);
    end
    for j = ch.H2   % a posteriori computed H2 norms
        A = Info.CL.a;
        B = Info.CL.b*ch.In{j};
        C = ch.Out{j}*Info.CL.c; 
        D = ch.Out{j}*Info.CL.d*ch.In{j};
        sys.CLch{j} = ss(A,B,C,D,opts.Ts);
        Info.H2(j) = norm(sys.CLch{j});        
        if opts.Ts == 0 % compensate for numerical errors in feedthrough
            if opts.Dc ~= 0
                if norm(DDzu{j}) > 0 
                    if norm(DDyw{j}) > 0 
                        sys.CLch{j}.d = zeros(nnz{j},nnw{j});
                        Info.H2(j) = norm(sys.CLch{j});
                    end
                end
            end
        end
    end
    if opts.Ts == 0 % check stability of closed-loop A-matrix
        Info.stab = max(real(eig(Info.CL.a)));
    else
        Info.stab = max(abs(eig(Info.CL.a)));
    end
end
 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Info] = Gshaping(A,Bw,Bu,Cz,Cy,Dzw,Dzu,Dyw,Dyu,a,b,ch,opts)

tol     = 1e-4; % tolerance primal residual
tolsing = 1e-9; % tolerance algebraic loop

% dimensions
nx = size(A,1); nu = size(Bu,2); ny = size(Cy,1);

% generalized plant
sys.P = ss(A,[Bw,Bu],[Cz;Cy],[Dzw,Dzu;Dyw,Dyu],opts.Ts);

% Total number of channels
nocH2   = length(ch.H2);
nocHinf = length(ch.Hinf);
noc     = nocH2 + nocHinf;

% System matrices and dimensions for specific channels
[BBw,CCz,DDzw,DDzu,DDyw,nnw,nnz] = deal(cell(1,noc));
for j=1:noc
    BBw{j}  = Bw*ch.In{j};
    CCz{j}  = ch.Out{j}*Cz;
    DDzw{j} = ch.Out{j}*Dzw*ch.In{j};
    DDzu{j} = ch.Out{j}*Dzu;
    DDyw{j} = Dyw*ch.In{j};
    nnw{j}  = size(BBw{j},2);
    nnz{j}  = size(CCz{j},1);
end

% new LMI system
LMIs = [];

% LMI variables
X = sdpvar(nx,nx,'full'); % slack matrix variables
Y = sdpvar(nx,nx,'full');
S = sdpvar(nx,nx,'full');
Q = sdpvar(nx,nx,'full'); % transformed controller variables
F = sdpvar(nx,ny,'full');
L = sdpvar(nu,nx,'full');
R = sdpvar(nu,ny,'full');
if opts.Dc == 0;          % strictly proper controller
    R = zeros(nu,ny);
end
[P,J,H] = deal(cell(1,noc)); % Lyapunov matrices for all channels
for j = 1:noc 
    P{j} = sdpvar(nx,nx,'symmetric');
    J{j} = sdpvar(nx,nx,'full');
    H{j} = sdpvar(nx,nx,'symmetric');
end
gam = cell(1,length(ch.Hinf));
for j = ch.Hinf
    gam{j} = sdpvar(1); % Hinf bounds
end
[W,mu] = deal(cell(1,length(ch.H2)));
for j = ch.H2
    W{j}  = sdpvar(nnz{j},nnz{j},'symmetric');
    mu{j} = sdpvar(1);  % H2 bounds
end

Info.LMIrows = 0;
for j = ch.Hinf % H-infinity LMIs
    T11 = P{j};
    T12 = J{j};
    T13 = A*X+Bu*L;
    T14 = A+Bu*R*Cy;
    T15 = BBw{j}+Bu*R*DDyw{j};
    T16 = zeros(nx,nnz{j});
    T22 = H{j};
    T23 = Q;
    T24 = Y*A+F*Cy;
    T25 = Y*BBw{j}+F*DDyw{j};
    T26 = zeros(nx,nnz{j});
    T33 = X+X'-P{j};
    T34 = eye(nx)+S'-J{j};
    T35 = zeros(nx,nnw{j});
    T36 = (CCz{j}*X+DDzu{j}*L)';
    T44 = Y+Y'-H{j};
    T45 = zeros(nx,nnw{j});
    T46 = (CCz{j}+DDzu{j}*R*Cy)';
    T55 = gam{j}*eye(nnw{j});
    T56 = (DDzw{j}+DDzu{j}*R*DDyw{j})';
    T66 = gam{j}*eye(nnz{j});
    
    Term = [T11 , T12 , T13 , T14 , T15 , T16;
        T12', T22 , T23 , T24 , T25 , T26;
        T13', T23', T33 , T34 , T35 , T36;
        T14', T24', T34', T44 , T45 , T46;
        T15', T25', T35', T45', T55 , T56;
        T16', T26', T36', T46', T56', T66];
    
    LMIs = [LMIs, Term >= 0];
    Info.LMIrows = Info.LMIrows + 4*nx + nnw{j} + nnz{j};
end

for j = ch.H2 % H2 LMIs
    T11 = P{j};
    T12 = J{j};
    T13 = A*X+Bu*L;
    T14 = A+Bu*R*Cy;
    T15 = BBw{j}+Bu*R*DDyw{j};
    T22 = H{j};
    T23 = Q;
    T24 = Y*A+F*Cy;
    T25 = Y*BBw{j}+F*DDyw{j};
    T33 = X+X'-P{j};
    T34 = eye(nx)+S'-J{j};
    T35 = zeros(nx,nnw{j});
    T44 = Y+Y'-H{j};
    T45 = zeros(nx,nnw{j});
    T55 = mu{j}*eye(nnw{j});
    
    Term = [T11 , T12 , T13 , T14 , T15;
            T12', T22 , T23 , T24 , T25;
            T13', T23', T33 , T34 , T35;
            T14', T24', T34', T44 , T45;
            T15', T25', T35', T45', T55];
    
    LMIs = [LMIs, Term >= 0];
    Info.LMIrows = Info.LMIrows + 4*nx + nnw{j};
    
    T11 = W{j};
    T12 = CCz{j}*X+DDzu{j}*L;
    T13 = CCz{j}+DDzu{j}*R*Cy;
    T14 = DDzw{j}+DDzu{j}*R*DDyw{j};
    T22 = X+X'-P{j};
    T23 = eye(nx)+S'-J{j};
    T24 = zeros(nx,nnw{j});
    T33 = Y+Y'-H{j};
    T34 = zeros(nx,nnw{j});
    T44 = mu{j}*eye(nnw{j});
    
    Term = [T11 , T12 , T13 , T14;
            T12', T22 , T23 , T24;
            T13', T23', T33 , T34;
            T14', T24', T34', T44];
    
    LMIs = [LMIs, Term >= 0];
    Info.LMIrows = Info.LMIrows + nnz{j} + 2*nx + nnw{j};
    
    LMIs = [LMIs, trace(W{j}) <= mu{j}];
end

% constraints
for j = ch.Hinf 
    if b(j) > 0 % Hinf bound channel j
        LMIs = [LMIs, gam{j} <= b(j)];
    end
end
for j = ch.H2 
    if b(j) > 0 % H2 bound channel j
        LMIs = [LMIs, mu{j} <= b(j)];
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
sol = solvesdp(LMIs,obj,sdpsettings('verbose',opts.verbose,'solver',opts.solver));
Info.cpusec = sol.solvertime;
p_res = checkset(LMIs);

% extract solution when it exists
Info.feas = 0;
if p_res > -tol
    Info.feas = 1;
    
    % solution variables
    X = double(X); 
    Y = double(Y);
    S = double(S);
    Q = double(Q); 
    F = double(F); 
    L = double(L); 
    if opts.Dc ~= 0
        R = double(R);
    end
    for j = 1:noc
        P{j} = double(P{j});
        J{j} = double(J{j});
        H{j} = double(H{j});
        Info.P{j} = [P{j}, J{j}; J{j}', H{j}];
    end
    
    % back-transformation of controller variables
    [Ac,Bc,Cc,Dc] = controller_transf(A,Bu,Cy,X,Y,Q,F,L,R,S); 
  
    % compensation for Dyu ~= 0
    if norm(Dyu,1) > 0
        if norm(Dc,1) > 0
            Myuc = eye(ny)+Dyu*Dc; 
            Mcyu = eye(nu)+Dc*Dyu;
            SV = svd(Myuc);
            if min(SV) < sqrt(tolsing)
                error('Algebraic loop');
            else
                Cc = Mcyu\Cc;
                Dc = Mcyu\Dc;
                Ac = Ac-Bc*Dyu*Cc;
                Bc = Bc/Myuc;
            end
        else
            Ac = Ac-Bc*Dyu*Cc;
        end
    end
    Info.K = ss(Ac,Bc,Cc,Dc,opts.Ts);  % controller
    Info.CL = lft(sys.P,Info.K,nu,ny); % closed-loop system
    for j = ch.Hinf % Hinf upper bounds
        Info.gam(j) = double(gam{j});
    end
    for j = ch.H2 % H2 upper bounds
        Info.mu(j) = double(mu{j});
    end
    sys.CLch = cell(1,noc);
    for j = ch.Hinf % a posteriori computed Hinf norms
        A = Info.CL.a;
        B = Info.CL.b*ch.In{j};
        C = ch.Out{j}*Info.CL.c; 
        D = ch.Out{j}*Info.CL.d*ch.In{j};
        sys.CLch{j} = ss(A,B,C,D,opts.Ts);
        Info.Hinf(j) = norm(sys.CLch{j},inf);
    end
    for j = ch.H2 % a posteriori computed H2 norms
        A = Info.CL.a;
        B = Info.CL.b*ch.In{j};
        C = ch.Out{j}*Info.CL.c; 
        D = ch.Out{j}*Info.CL.d*ch.In{j};
        sys.CLch{j} = ss(A,B,C,D,opts.Ts);
        Info.H2(j) = norm(sys.CLch{j});
    end   
    if opts.Ts == 0 % check stability of closed-loop A-matrix
        Info.stab = max(real(eig(Info.CL.a)));
    else
        Info.stab = max(abs(eig(Info.CL.a)));
    end    
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Ac,Bc,Cc,Dc] = controller_transf(A,Bu,Cy,X,Y,Q,F,L,R,S)
% Dimensions
nx = size(A,1); nu = size(Bu,2); ny = size(Cy,1);                    

% Apply singular value decomposition
[Left,D1,Right] = svd(S-Y*X); 
V = Left*sqrt(D1); U = sqrt(D1)*Right';

% Transformation
K = [ V^(-1), -V^(-1)*Y*Bu; zeros(nu,nx), eye(nu)]*...
    [Q-Y*A*X, F; L, R]*...
    [U^(-1), zeros(nx,ny); -Cy*X*U^(-1), eye(ny)];

% Extract controller state-space matrices
Ac = K(1:nx    , 1:nx); Bc = K(1:nx    , nx+1:end);
Cc = K(nx+1:end, 1:nx); Dc = K(nx+1:end, nx+1:end);
end
