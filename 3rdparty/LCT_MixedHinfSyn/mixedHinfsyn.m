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


% Controller design subject to multiple Hinf specifications on closed-loop
% transfer functions from a common exogenous input w to various exogenous
% outputs z_i.
%
% Inputs
%   sys_P : generalized plant with inputs [w;u] and outputs [z_i;y] 
%   my    : dimension of the measured output y
%   mu    : dimension of the control input u
%   Mz    : dimensions of regulated outputs z_i: [mz_i]
%   alpha : weigths determining how to treat the Hinf constraints
%           alpha_i = 0 -> ||H_i||_inf < 1
%           alpha_i > 0 -> ||H_i||_inf < gamma_i
%           objective = sum(alpha_i*gamma_i^2)
%
% Outputs
%   sys_K : optimal controller
%   sys_H : corresponding closed-loop system 
%   gamma : Hinf norms of sys_H: [gamma_i_opt, gamma_i_actual]
%
% Goele Pipeleers, 2009
%
%
% References
%   P. Gahinet, Explicit controller formulas for LMI-based H_\infty 
%   synthesis. Automatica, 32(7):1007-1014, 1996.
%   P. Gahinet and P. Apkarian, A linear matrix inequality appoach to
%   H_\infty control. International Journal of Robust and Nonlinear
%   Control, 4:421-448, 1994.
%   T. Iwasaki and R.E. Skelton, All controllers for te general H_\infty
%   control problem: LMI existence conditions and state-space formulas.
%   Automatica, 30(8):1307-1317, 1994.



function [sys_K, sys_H, gamma] = mixedHinfsyn (sys_P, my, mu, Mz, alpha, options)


% STEP 0. Default options
% -----------------------
if nargin < 5
    error('Input arguments incomplete')
elseif nargin == 5
    options.gamma.solver = 'lmilab';
    options.gamma.solution = 2;
    options.controller.solution = 1;
    options.controller.LMItype = 1;
    options.controller.LMIsol = 2;
end


% STEP 1. Compute the (sub)optimal gamma_i
% ----------------------------------------

[gamma, X, Y, Gamma] = mixedHinfsyn_gamma (sys_P, my, mu, Mz, alpha, options);



% STEP 2. Compute the corresponding controller
% --------------------------------------------
%     => currently only full order!!

if sys_P.Ts == 0    % -> cont.
    if options.controller.solution == 1
        [sys_K, sys_H, check_LMI] = mixedHinfsyn_K1c(sys_P, my, mu, X, Y, Gamma, options);
    elseif options.controller.solution == 2
        [sys_K, sys_H, check_LMI] = mixedHinfsyn_K2c(sys_P, my, mu, X, Y, Gamma, options);
    else
        error('unvalid entry for options.controller.solution')
    end
    
else                % -> discr. 
    if options.controller.solution == 1
        [sys_K, sys_H, check_LMI] = mixedHinfsyn_K1d(sys_P, my, mu, X, Y, Gamma, options);
    elseif options.controller.solution == 2
        [sys_K, sys_H, check_LMI] = mixedHinfsyn_K2d(sys_P, my, mu, X, Y, Gamma, options);
    else
        error('unvalid entry for options.controller.solution')
    end
end



% STEP 3. Analysis of closed-loop system
% --------------------------------------

noc = length(alpha);
Iz2 = cumsum(Mz);
Iz1 = [1;Iz2(1:end-1)+1];

gamma_act = zeros(noc,1);
for i = 1:noc
   gamma_act(i) =  norm(sys_H(Iz1(i):Iz2(i),:), inf);
end

gamma = [gamma, gamma_act];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [gamma, X, Y, Gamma] = mixedHinfsyn_gamma (sys_P, my, mu, Mz, alpha, options)


% STEP 0. Settings
% ----------------

if strcmp(options.gamma.solver, 'lmilab')
    LMIparser = 2;                      % 2: LMI Lab
else
    LMIparser = 1;                      % 1: yalmip
end

solution = options.gamma.solution;      % 1: scaled Hinf norm
                                        % 2: Lyap. shaping



% STEP 1. State-space model of generalized plant
% ----------------------------------------------

% 1.0. Continuous or discrete time
if sys_P.Ts == 0    % -> cont.
    Phi = [0,1;1,0];
else                % -> discr. 
    Phi = [1,0;0,-1];
end


% 1.1. State-space model and dimensions
[A,B,C,D] = ssdata(sys_P);
n = length(A(:,1));
mz = sum(Mz);
mw = length(B(1,:))-mu;


% 1.2. Subdivision according to inputs and outputs
                        Bw = B(:,1:mw);             Bu = B(:,mw+1:end);
Cz = C(1:mz,:);         Dzw = D(1:mz,1:mw);         Dzu = D(1:mz,mw+1:end);
Cy = C(mz+1:end,:);     Dyw = D(mz+1:end,1:mw);


% 1.3. Various Hinf constraints
if length(Mz) ~= length(alpha)
    error('Dimensions of Mz and alpha are inconsistent')
end
if sum(alpha<0)~=0
    error('weights must be positive')
end
if my+mz ~= length(C(:,1))
    error('Dimensions of y and z are inconsistent')
end
Mz = Mz(:);
alpha = alpha(:);
noc = length(alpha);
Iz2 = cumsum(Mz);
Iz1 = [1;Iz2(1:end-1)+1];



% STEP 2. SDP to compute the optimal gamma_i
% ------------------------------------------

% 2.0. Preprocessing: null-spaces
V = null([Cy, Dyw]);
[rV,cV] = size(V);
V = [V, zeros(rV,mz); zeros(mz,cV), eye(mz)];

W = null([Bu', Dzu']);
[rW,cW] = size(W);
W = [W, zeros(rW,mw); zeros(mw,cW), eye(mw)];


% 2.1. Parsing the SDP with Yalmip    
switch LMIparser
    case 1  % -> yalmip

%   2.1.1. Declaring the variables        
X = sdpvar(n);      % nxn top-left subblock of P
Y = sdpvar(n);      % nxn top-left subblock of inv(P)
switch solution
    case 1  % -> scaled Hinf
        Gamma = diag(sdpvar(mz,1));
    case 2  % -> Lyap. shaping
        Gamma = sdpvar(mz);
end
c = zeros(mz,1);
for i = 1:noc
    c(Iz1(i):Iz2(i)) = alpha(i)/Mz(i);
    if alpha(i) == 0
        Gamma(Iz1(i):Iz2(i),Iz1(i):Iz2(i)) = eye(Mz(i));
    else
        Gamma(Iz1(i):Iz2(i),Iz1(i):Iz2(i)) = sdpvar(1) * eye(Mz(i));
    end
end

%   2.1.2. Specifying the LMIs
ABw = [A, Bw, zeros(n,mz); eye(n), zeros(n,mw+mz)];
Zx = ABw'*kron(Phi,X)*ABw + [zeros(n,n+mw),Cz'; zeros(mw,n),-eye(mw),Dzw'; Cz,Dzw,-Gamma];
ACz = [A', Cz', zeros(n,mw); eye(n), zeros(n,mz+mw)];
Zy = ACz'*kron(Phi,Y)*ACz + [zeros(n,n+mz),Bw; zeros(mz,n),-Gamma,Dzw; Bw',Dzw',-eye(mw)];
constr = ([X,eye(n);eye(n),Y] >= 0) + (V'*Zx*V <= 0) + (W'*Zy*W <= 0);

%   2.1.3. Specifying the objective
goal = trace(diag(c)*Gamma);

%   2.1.4. Solving the SDP
options = sdpsettings('solver',options.gamma.solver);
sol = solvesdp(constr, goal, options);
check_LMIs = {'P pos' , min(eig(double([X,eye(n);eye(n),Y])));
              'Zx neg', max(eig(double(V'*Zx*V)))            ;
              'Zy neg', max(eig(double(W'*Zy*W)))            }
X = double(X);
Y = double(Y);
Gamma = double(Gamma);
gamma = diag(Gamma);
gamma = sqrt(gamma(Iz1));


% 2.2. Parsing the SDP with LMI Lab
	case 2  % -> LMI Lab

%   2.2.1. Declaring the variables
setlmis([]);
X = lmivar(1, [n,1]);       % nxn top-left subblock of P
[Y,mx] = lmivar(1, [n,1]);	% nxn top-left subblock of inv(P)
Gamma0 = zeros(mz);
strucGamma = zeros(mz);
for i = 1:noc
    if alpha(i) == 0
        Gamma0(Iz1(i):Iz2(i),Iz1(i):Iz2(i)) = eye(Mz(i));
    else
        strucGamma(Iz1(i):Iz2(i),Iz1(i):Iz2(i)) = ((mx+1)/2) *eye(Mz(i));
        mx = mx+1;
    end
    if solution == 2
        temp = reshape(mx+(1:Mz(i)*sum(Mz(i+1:end)))', Mz(i), sum(Mz(i+1:end)));
        strucGamma(Iz1(i):Iz2(i),Iz2(i)+1:end) = temp;
        mx = mx + Mz(i)*sum(Mz(i+1:end));
    end
end
strucGamma = strucGamma+strucGamma';
Gamma = lmivar(3, strucGamma);


%   2.2.2. Specifying the LMIs
ABw = [A, Bw];
IO = [eye(n), zeros(n,mw)];
OI = [zeros(mw,n), eye(mw)];
Zx_neg = newlmi;
lmiterm([Zx_neg,1,1,X], Phi(1,1)*ABw', ABw);
lmiterm([Zx_neg,1,1,X], Phi(1,2)*ABw', IO, 's');
lmiterm([Zx_neg,1,1,X], Phi(2,2)*IO', IO);
lmiterm([Zx_neg,1,1,0], -OI'*eye(mw)*OI);
lmiterm([Zx_neg,1,2,0], [Cz';Dzw']);
if ~isempty(Gamma)
    lmiterm([Zx_neg,2,2,Gamma], -1, 1);
end
lmiterm([Zx_neg,2,2,0], -Gamma0);
lmiterm([Zx_neg,0,0,0],V);

ACz = [A', Cz'];
IO = [eye(n), zeros(n,mz)];
OI = [zeros(mz,n), eye(mz)];
Zy_neg = newlmi;
lmiterm([Zy_neg,1,1,Y], Phi(1,1)*ACz', ACz);
lmiterm([Zy_neg,1,1,Y], Phi(1,2)*ACz', IO, 's');
lmiterm([Zy_neg,1,1,Y], Phi(2,2)*IO', IO);
if ~isempty(Gamma)
    lmiterm([Zy_neg,1,1,Gamma], -OI', OI);
end
lmiterm([Zy_neg,1,1,0], -OI'*Gamma0*OI);
lmiterm([Zy_neg,1,2,0], [Bw;Dzw]);
lmiterm([Zy_neg,2,2,0], -eye(mw));
lmiterm([Zy_neg,0,0,0],W);

P_pos = newlmi;
lmiterm([-P_pos,1,1,X], 1, 1);
lmiterm([-P_pos,1,2,0], 1);
lmiterm([-P_pos,2,2,Y], 1, 1);

LMIs = getlmis;

%   2.2.3. Specifying the objective
c = zeros(decnbr(LMIs),1);
if ~isempty(Gamma)
    diagGamma = diag(decinfo(LMIs,Gamma));
    diagGamma = diagGamma(Iz1);
    c(diagGamma(diagGamma~=0)) = alpha((diagGamma~=0));
end
    
%   2.2.4. Solving the SDP
[fopt, xopt] = mincx(LMIs, c);
LMIsopt = evallmi(LMIs, xopt);
[Zx_neg, trash] = showlmi(LMIsopt, Zx_neg);
[Zy_neg, trash] = showlmi(LMIsopt, Zy_neg);
[trash, P_pos] = showlmi(LMIsopt, P_pos);
check_LMIs = {'P pos' , min(eig(P_pos)) ;
              'Zx neg', max(eig(Zx_neg));
              'Zy neg', max(eig(Zy_neg))}

X = dec2mat(LMIs, xopt, X);
Y = dec2mat(LMIs, xopt, Y);
Gamma = Gamma0 + dec2mat(LMIs, xopt, Gamma);
gamma = diag(double(Gamma));
gamma = sqrt(gamma(Iz1));

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [sys_K, sys_H, check_LMI] = mixedHinfsyn_K1c(sys_P, my, mu, X, Y, Gamma, options)


% STEP 0. Settings
% ----------------

LMItype = options.controller.LMItype;       % 1: LMI in terms of transformed controller variables
                                            % 2: LMI in terms of actual controller variables
                                            %       (!?! usually ill-conditioned -> solve with lmilab or analytical)    
LMIsol = options.controller.LMIsol;         % 1: yalmip
                                            % 2: 'basiclmi' from robust control toolbox
tolsing = 1e-9;                   



% STEP 1. State-space model of generalized plant
% ----------------------------------------------

% 1.1. State-space model and dimensions
[A,B,C,D] = ssdata(sys_P);
n = length(A(:,1));
mw = length(B(1,:))-mu;
mz = length(C(:,1))-my;


% 1.2. Subdivision according to inputs and outputs
                        Bw = B(:,1:mw);             Bu = B(:,mw+1:end);
Cz = C(1:mz,:);         Dzw = D(1:mz,1:mw);         Dzu = D(1:mz,mw+1:end);
Cy = C(mz+1:end,:);     Dyw = D(mz+1:end,1:mw);     Dyu = D(mz+1:end,mw+1:end);



% STEP 2. Compute M and N by singular value decomposition
% -------------------------------------------------------

[U,S,V] = svd(eye(n)-Y*X);
N = U*sqrt(S);
M = V*sqrt(S);



% STEP 3. Solve LMI feasibility problem
% -------------------------------------

% 3.1. Formulate LMI: Z+U'*theta*V+V'*theta'*U < 0
switch LMItype
    case 1  % -> transformed controller parameters
Z = [A*Y+Y*A', A       , Bw      , Y*Cz' ; 
     A'      , X*A+A'*X, X*Bw    , Cz'   ;
     Bw'     , Bw'*X   , -eye(mw), Dzw'  ;
     Cz*Y    , Cz      , Dzw     , -Gamma];
U = [zeros(n), eye(n)     , zeros(n,mw) , zeros(n,mz);
     Bu'     , zeros(mu,n), zeros(mu,mw), Dzu'       ];
V = [eye(n)     , zeros(n), zeros(n,mw), zeros(n,mz) ;
     zeros(my,n), Cy      , Dyw        , zeros(my,mz)];

    case 2  % -> actual controller parameters (!?! usually ill-conditioned -> solve with lmilab or analytical)
XAY = 0.5*(X*A*Y+(Y*A'*X)');
Z = [A*Y+Y*A', A+XAY'  , Bw      , Y*Cz' ; 
     A'+XAY  , X*A+A'*X, X*Bw    , Cz'   ;
     Bw'     , Bw'*X   , -eye(mw), Dzw'  ;
     Cz*Y    , Cz      , Dzw     , -Gamma];
U = [zeros(n), M'   , zeros(n,mw) , zeros(n,mz);
     Bu'     , Bu'*X, zeros(mu,mw), Dzu'       ];
V = [N'  , zeros(n), zeros(n,mw), zeros(n,mz) ;
     Cy*Y, Cy      , Dyw        , zeros(my,mz)];
end


% 3.2. Solve LMI feasibility problem
switch LMIsol
    case 1  % -> yalmip
        theta = sdpvar(n+mu,n+my,'full');
        eps = sdpvar(1);
        constr = (Z+U'*theta*V+V'*theta'*U <= eps*eye(2*n+mw+mz));
        goal = eps;
        options = sdpsettings('solver','sdpt3');
        sol = solvesdp(constr, goal, options);
        theta = double(theta);
        double(eps)

    case 2  % -> 'basiclmi' from robust control toolbox
        theta = basiclmi(Z,U,V,'Xmin,Shift');
end
check_LMI = max(eig(Z+U'*theta*V+V'*theta'*U))

 
% 3.3. Invert bijecive transformation of controller parameters
if LMItype == 1
    theta = theta - [X*A*Y, zeros(n,my); zeros(mu,n+my)];
    theta = [M, X*Bu; zeros(mu,n), eye(mu)] \ theta;
    theta = theta / [N', zeros(n,my); Cy*Y, eye(my)];
end

Ak = theta(1:n,1:n);
Bk = theta(1:n,n+1:end);
Ck = theta(n+1:end,1:n);
Dk = theta(n+1:end,n+1:end);



% STEP 4. Postprocessing
% ----------------------

% 4.1. Compensation for Dyu~=0
Ak0 = Ak;   Bk0 = Bk;   Ck0 = Ck;   Dk0 = Dk;

if norm(Dyu,1) > 0,
    if norm(Dk,1) > 0,
        Myuk = eye(my)+Dyu*Dk; 
        Mkyu = eye(mu)+Dk*Dyu;
        S = svd(Myuk);
        if min(S) < sqrt(tolsing)
            error('Algebraic loop');
        else
            Ck = Mkyu\Ck;
            Dk = Mkyu\Dk;
            Ak = Ak-Bk*Dyu*Ck;
            Bk = Bk/Myuk;
        end
    else
        Ak = Ak-Bk*Dyu*Ck;
    end
end
sys_K = ss(Ak,Bk,Ck,Dk);


% 4.2. Corresponding closed-loop system
sys_H = lft(sys_P, sys_K, mu, my);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [sys_K, sys_H, check_LMI] = mixedHinfsyn_K1d(sys_P, my, mu, X, Y, Gamma, options)


% STEP 0. Settings
% ----------------

LMItype = options.controller.LMItype;       % 1: LMI in terms of transformed controller variables
                                            % 2: LMI in terms of actual controller variables
                                            %       (!?! usually ill-conditioned -> solve with lmilab or analytical)
LMIsol = options.controller.LMIsol;         % 1: yalmip
                                            % 2: 'basiclmi' from robust
                                            % control toolbox
tolsing = 1e-9;                   



% STEP 1. State-space model of generalized plant
% ----------------------------------------------

% 1.1. State-space model and dimensions
[A,B,C,D] = ssdata(sys_P);
n = length(A(:,1));
mw = length(B(1,:))-mu;
mz = length(C(:,1))-my;


% 1.2. Subdivision according to inputs and outputs
                        Bw = B(:,1:mw);             Bu = B(:,mw+1:end);
Cz = C(1:mz,:);         Dzw = D(1:mz,1:mw);         Dzu = D(1:mz,mw+1:end);
Cy = C(mz+1:end,:);     Dyw = D(mz+1:end,1:mw);     Dyu = D(mz+1:end,mw+1:end);



% STEP 2. Compute M and N by singular value decomposition
% -------------------------------------------------------

[U,S,V] = svd(eye(n)-Y*X);
N = U*sqrt(S);
M = V*sqrt(S);



% STEP 3. Solve LMI feasibility problem
% -------------------------------------

% 3.1. Formulate LMI: Z+U'*theta*V+V'*theta'*U < 0
switch LMItype
    case 1  % -> transformed controller parameters
Z = [-Y         , -eye(n)    , A*Y        , A          , Bw         , zeros(n,mz);
     -eye(n)    , -X         , zeros(n)   , X*A        , X*Bw       , zeros(n,mz);
     Y*A'       , zeros(n)   , -Y         ,-eye(n)     , zeros(n,mw), Y*Cz'      ;
     A'         , A'*X       , -eye(n)    , -X         , zeros(n,mw), Cz'        ;
     Bw'        , Bw'*X      , zeros(mw,n), zeros(mw,n), -eye(mw)   , Dzw'       ;
     zeros(mz,n), zeros(mz,n), Cz*Y       , Cz         , Dzw        , -Gamma     ];   
U = [zeros(n), eye(n)     , zeros(n,2*n) , zeros(n,mw) , zeros(n,mz);
     Bu'     , zeros(mu,n), zeros(mu,2*n), zeros(mu,mw), Dzu'       ];
V = [zeros(n,2*n) , eye(n)     , zeros(n), zeros(n,mw), zeros(n,mz) ;
     zeros(my,2*n), zeros(my,n), Cy      , Dyw        , zeros(my,mz)];

    case 2  % -> actual controller parameters (!?! usually ill-conditioned -> solve with lmilab or analytical)
XAY = 0.5*(X*A*Y+(Y*A'*X)');
Z = [-Y         , -eye(n)    , A*Y        , A          , Bw         , zeros(n,mz);
     -eye(n)    , -X         , XAY        , X*A        , X*Bw       , zeros(n,mz);
     Y*A'       , XAY'       , -Y         ,-eye(n)     , zeros(n,mw), Y*Cz'      ;
     A'         , A'*X       , -eye(n)    , -X         , zeros(n,mw), Cz'        ;
     Bw'        , Bw'*X      , zeros(mw,n), zeros(mw,n), -eye(mw)   , Dzw'       ;
     zeros(mz,n), zeros(mz,n), Cz*Y       , Cz         , Dzw        , -Gamma     ];   
U = [zeros(n), M'   , zeros(n,2*n) , zeros(n,mw) , zeros(n,mz);
     Bu'     , Bu'*X, zeros(mu,2*n), zeros(mu,mw), Dzu'       ];
V = [zeros(n,2*n) , N'  , zeros(n), zeros(n,mw), zeros(n,mz) ;
     zeros(my,2*n), Cy*Y, Cy      , Dyw        , zeros(my,mz)];
end



% 3.2. Solve LMI feasibility problem
switch LMIsol
    case 1  % -> yalmip
        theta = sdpvar(n+mu,n+my,'full');
        eps = sdpvar(1);
        constr = (Z+U'*theta*V+V'*theta'*U <= eps*eye(4*n+mw+mz));
        goal = eps;
        options = sdpsettings('solver','sdpt3');
        sol = solvesdp(constr, goal, options);
        theta = double(theta);

    case 2  % -> 'basiclmi' from robust control toolbox
        theta = basiclmi(Z,U,V,'Xmin,Shift');
end
check_LMI = max(eig(Z+U'*theta*V+V'*theta'*U))

 
% 3.3. Invert bijecive transformation of controller parameters
if LMItype == 1
    theta = theta - [X*A*Y, zeros(n,my); zeros(mu,n+my)];
    theta = [M, X*Bu; zeros(mu,n), eye(mu)] \ theta;
    theta = theta / [N', zeros(n,my); Cy*Y, eye(my)];
end

Ak = theta(1:n,1:n);
Bk = theta(1:n,n+1:end);
Ck = theta(n+1:end,1:n);
Dk = theta(n+1:end,n+1:end);



% STEP 4. Postprocessing
% ----------------------

% 4.1. Compensation for Dyu~=0
Ak0 = Ak;   Bk0 = Bk;   Ck0 = Ck;   Dk0 = Dk;

if norm(Dyu,1) > 0,
    if norm(Dk,1) > 0,
        Myuk = eye(my)+Dyu*Dk; 
        Mkyu = eye(mu)+Dk*Dyu;
        S = svd(Myuk);
        if min(S) < sqrt(tolsing)
            error('Algebraic loop');
        else
            Ck = Mkyu\Ck;
            Dk = Mkyu\Dk;
            Ak = Ak-Bk*Dyu*Ck;
            Bk = Bk/Myuk;
        end
    else
        Ak = Ak-Bk*Dyu*Ck;
    end
end
sys_K = ss(Ak,Bk,Ck,Dk,sys_P.Ts);


% 4.2. Corresponding closed-loop system
sys_H = lft(sys_P, sys_K, mu, my);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [sys_K, sys_H, check_LMI] = mixedHinfsyn_K2c(sys_P, my, mu, X, Y, Gamma, options)


% STEP 0. Settings
% ----------------

LMIsol = options.controller.LMIsol;     % 1: yalmip
                                        % 2: 'basiclmi' from robust control toolbox
tolsing = 1e-9;                



% STEP 1. State-space model of generalized plant
% ----------------------------------------------

% 1.1. State-space model and dimensions
[A,B,C,D] = ssdata(sys_P);
n = length(A(:,1));
mw = length(B(1,:))-mu;
mz = length(C(:,1))-my;


% 1.2. Subdivision according to inputs and outputs
                        Bw = B(:,1:mw);             Bu = B(:,mw+1:end);
Cz = C(1:mz,:);         Dzw = D(1:mz,1:mw);         Dzu = D(1:mz,mw+1:end);
Cy = C(mz+1:end,:);     Dyw = D(mz+1:end,1:mw);     Dyu = D(mz+1:end,mw+1:end);



% STEP 2. Compute an adequate DDk
% -------------------------------
%   [-I, Dcl'; Dcl, -Gamma] < 0

% 2.1. Formulate LMI condition as: Z + U'*DDk*V + V'*DDk'*U < 0
Z = [-eye(mw), Dzw'; Dzw, -Gamma];
U = [zeros(mw,mu); Dzu]';
V = [Dyw, zeros(my,mz)];


% 2.2. Compute a feasible DDk
if max(eig(Z)) <= -1e-2;
    DDk = zeros(mu,my);
else
    switch LMIsol
        case 1  % -> yalmip
            DDk = sdpvar(mu,my,'full');
            eps = sdpvar(1);
            constr = (Z+U'*DDk*V+V'*DDk'*U <= eps*eye(mw+mz));
            goal = eps;
            options = sdpsettings('solver','sdpt3');
            sol = solvesdp(constr, goal, options);
            DDk = double(DDk);

        case 2  % -> 'basiclmi' from robust control toolbox
            DDk = basiclmi(Z,U,V,'Xmin,Shift');
    end
    if max(eig(Z)) <= max(eig(Z+U'*DDk*V+V'*DDk'*U))
        DDk = zeros(mu,my);
    end
end
check_LMI_DDk = max(eig(Z+U'*DDk*V+V'*DDk'*U))


% 2.3. Post-processing: auxiliary matrices
AA = A + Bu*DDk*Cy;
BBw = Bw + Bu*DDk*Dyw;
CCz = Cz + Dzu*DDk*Cy;
Dcl = Dzw + Dzu*DDk*Dyw;

Delta = -(Z+U'*DDk*V+V'*DDk'*U);
cholDelta = chol(Delta)';           % Delta = cholDelta * cholDelta'



% STEP 3. Compute BBk
% -------------------

% 3.1. Compute regular contribution
if rank(Dyw,tolsing) == 0
    BBk = zeros(n,my);
else
    Z = [zeros(my)   , Dyw    , zeros(my,mz);
         Dyw'        ,-eye(mw), Dcl'        ;
         zeros(mz,my), Dcl    ,-Gamma       ];
    U = - [Cy; Bw'*X; CCz];
    BBk = Z\U;  % pinv(Z)*U;
    BBk = BBk(1:my,:)';
end


% 3.2. Detect singularity
[U,S,V] = svd(Dyw);
S = S(1:min(my,mw), 1:min(my,mw));
i2 = find(diag(S) <= tolsing);
U2 = U(:,i2);
Pi_yw = U2*U2';
sing_yw = (norm(U2'*Cy) > tolsing*norm(Cy,1));


% 3.3. Additional contribution for singular case
if sing_yw
    BBk1 = BBk;

    aux1 = A'*X+Cy'*BBk1';
    aux2 = cholDelta \ [Bw'*X+Dyw'*BBk1'; CCz];
    Z = aux1+aux1' + aux2'*aux2;
    V = Pi_yw*Cy;
    U = eye(n);
    switch LMIsol
        case 1  % -> yalmip
            BBk2 = sdpvar(n,my,'full');
            eps = sdpvar(1);
            constr = (Z+U'*BBk2*V+V'*BBk2'*U <= eps*eye(n));
            goal = eps;
            options = sdpsettings('solver','sdpt3');
            sol = solvesdp(constr, goal, options);
            BBk2 = double(BBk2);

        case 2  % -> 'basiclmi' from robust control toolbox
            BBk2 = basiclmi(Z,U,V,'Xmin,Shift');
    end
    if max(eig(Z)) <= max(eig(Z+U'*BBk2*V+V'*BBk2'*U))
        BBk2 = zeros(n,my);
    end
    
    BBk = BBk1 + BBk2*Pi_yw;
end


% 3.4. Verification of the result
Z = [A'*X+X*A, X*Bw    , CCz'  ;
     Bw'*X   , -eye(mw), Dcl'  ;
     CCz     , Dcl     , -Gamma];
U = [eye(n), zeros(n,mw+mz)];
V = [Cy, Dyw, zeros(my,mz)];
check_LMI_BBk = max(eig(Z+U'*BBk*V+V'*BBk'*U))



% STEP 4. Compute CCk
% -------------------

% 4.1. Compute regular contribution
if rank(Dzu,tolsing) == 0
    CCk = zeros(mu,n);
else
    Z = [zeros(mu)   , zeros(mu,mw), Dzu' ;
         zeros(mw,mu),-eye(mw)     , Dcl' ;
         Dzu         , Dcl         ,-Gamma];
    U = - [Bu'; BBw'; Cz*Y];
    CCk = Z\U;
    CCk = CCk(1:mu,:);
end


% 4.2. Detect singularity: (I-pinv(Dzu)*Dzu)*Bu' = 0
[U,S,V] = svd(Dzu);
S = S(1:min(mz,mu), 1:min(mz,mu));
i2 = find(diag(S) <= tolsing);
V2 = V(:,i2);
Pi_zu = V2*V2';
sing_zu = (norm(V2'*Bu') > tolsing*norm(Bu,1))


% 4.3. Additional contribution for singular case
if sing_zu
    CCk1 = CCk;

    aux1 = A*Y+Bu*CCk1;
    aux2 = cholDelta \ [BBw'; Cz*Y+Dzu*CCk1];
    Z = aux1+aux1' + aux2'*aux2;
    U = (Bu*Pi_zu)';
    V = eye(n);
    switch LMIsol
        case 1  % -> yalmip
            CCk2 = sdpvar(mu,n,'full');
            eps = sdpvar(1);
            constr = (Z+U'*CCk2*V+V'*CCk2'*U <= eps*eye(n));
            goal = eps;
            options = sdpsettings('solver','sdpt3');
            sol = solvesdp(constr, goal, options);
            CCk2 = double(CCk2);

        case 2  % -> 'basiclmi' from robust control toolbox
            CCk2 = basiclmi(Z,U,V,'Xmin,Shift');
    end
    if max(eig(Z)) <= max(eig(Z+U'*CCk2*V+V'*CCk2'*U))
        CCk2 = zeros(n,my);
    end

    CCk = CCk1 + Pi_zu*CCk2;
end


% 4.4. Verification of te result
Z = [A*Y+Y*A', BBw     , Y*Cz' ;
     BBw'    , -eye(mw), Dcl'  ;
     Cz*Y    , Dcl     , -Gamma];
U = [Bu', zeros(mu,mw), Dzu'];
V = [eye(n), zeros(n,mw+mz)];
check_LMI_CCk = max(eig(Z+U'*CCk*V+V'*CCk'*U))



% STEP 5. Compute AAk
% -------------------

aux1 = cholDelta \ [Bw'*X+Dyw'*BBk'; CCz];
aux2 = cholDelta \ [BBw'; Cz*Y+Dzu*CCk];
AAk = - (AA' + aux1'*aux2);



% STEP 6. Invert bijecive transformation of controller parameters
% ---------------------------------------------------------------

% 6.1. Compute M and N by singular value decomposition
[U,S,V] = svd(eye(n)-Y*X);
N = U*sqrt(S);
M = V*sqrt(S);


% 6.2. Reconstructing the controller
theta = [AAk,BBk;CCk,DDk] - [X*A*Y, zeros(n,my); zeros(mu,n+my)];
theta = [M, X*Bu; zeros(mu,n), eye(mu)] \ theta;
theta = theta / [N', zeros(n,my); Cy*Y, eye(my)];

Ak = theta(1:n,1:n);
Bk = theta(1:n,n+1:end);
Ck = theta(n+1:end,1:n);
Dk = theta(n+1:end,n+1:end);



% STEP 7. Verification of the result
% ----------------------------------

% 7.1. LMI in terms of the transformed controller parameters
Z = [A*Y+Y*A', A       , Bw      , Y*Cz' ; 
     A'      , X*A+A'*X, X*Bw    , Cz'   ;
     Bw'     , Bw'*X   , -eye(mw), Dzw'  ;
     Cz*Y    , Cz      , Dzw     , -Gamma];
U = [zeros(n), eye(n)     , zeros(n,mw) , zeros(n,mz);
     Bu'     , zeros(mu,n), zeros(mu,mw), Dzu'       ];
V = [eye(n)     , zeros(n), zeros(n,mw), zeros(n,mz) ;
     zeros(my,n), Cy      , Dyw        , zeros(my,mz)];
theta = [AAk,BBk;CCk,DDk]
check_LMI1 = max(eig(Z+U'*theta*V+V'*theta'*U));


% 7.2. LMI in terms of the controller parameters
Z = [A*Y+Y*A', A+Y*A'*X, Bw      , Y*Cz' ; 
     A'+X*A*Y, X*A+A'*X, X*Bw    , Cz'   ;
     Bw'     , Bw'*X   , -eye(mw), Dzw'  ;
     Cz*Y    , Cz      , Dzw     , -Gamma];
U = [zeros(n), M'   , zeros(n,mw) , zeros(n,mz);
     Bu'     , Bu'*X, zeros(mu,mw), Dzu'       ];
V = [N'  , zeros(n), zeros(n,mw), zeros(n,mz) ;
     Cy*Y, Cy      , Dyw        , zeros(my,mz)];
theta = [Ak,Bk;Ck,Dk]; 
check_LMI2 = max(eig(Z+U'*theta*V+V'*theta'*U));

check_LMI = [check_LMI1;check_LMI2]



% STEP 8. Postprocessing
% ----------------------

% 8.1. Compensation for Dyu~=0
Ak0 = Ak;   Bk0 = Bk;   Ck0 = Ck;   Dk0 = Dk;

if norm(Dyu,1) > 0,
    if norm(Dk,1) > 0,
        Myuk = eye(my)+Dyu*Dk; 
        Mkyu = eye(mu)+Dk*Dyu;
        S = svd(Myuk);
        if min(S) < sqrt(tolsing)
            error('Algebraic loop');
        else
            Ck = Mkyu\Ck;
            Dk = Mkyu\Dk;
            Ak = Ak-Bk*Dyu*Ck;
            Bk = Bk/Myuk;
        end
    else
        Ak = Ak-Bk*Dyu*Ck;
    end
end
sys_K = ss(Ak,Bk,Ck,Dk);


% 8.2. Corresponding closed-loop system
sys_H = lft(sys_P, sys_K, mu, my);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [sys_K, sys_H, check_LMI] = mixedHinfsyn_K2d(sys_P, my, mu, X, Y, Gamma, options)


% STEP 0. Settings
% ----------------

LMIsol = options.controller.LMIsol;     % 1: yalmip
                                        % 2: 'basiclmi' from robust control toolbox
tolsing = 1e-9;                



% STEP 1. State-space model of generalized plant
% ----------------------------------------------

% 1.1. State-space model and dimensions
[A,B,C,D] = ssdata(sys_P);
n = length(A(:,1));
mw = length(B(1,:))-mu;
mz = length(C(:,1))-my;


% 1.2. Subdivision according to inputs and outputs
                        Bw = B(:,1:mw);             Bu = B(:,mw+1:end);
Cz = C(1:mz,:);         Dzw = D(1:mz,1:mw);         Dzu = D(1:mz,mw+1:end);
Cy = C(mz+1:end,:);     Dyw = D(mz+1:end,1:mw);     Dyu = D(mz+1:end,mw+1:end);



% STEP 2. Compute an adequate DDk
% -------------------------------

% 2.1. Formulate LMI condition as: Z + U'*DDk*V + V'*DDk'*U < 0
Z = [-Y         , A          , Bw         , zeros(n,mz);
     A'         , -X         , zeros(n,mw), Cz'        ;
     Bw'        , zeros(mw,n), -eye(mw)   , Dzw'       ;
     zeros(mz,n), Cz         , Dzw        , -Gamma     ];
U = [Bu', zeros(mu,n), zeros(mu,mw), Dzu'];
V = [zeros(my,n), Cy, Dyw, zeros(my,mz)];


% 2.2. Compute a feasible DDk
if max(eig(Z)) <= -1e-2;
    DDk = zeros(mu,my);
else
    switch LMIsol
        case 1  % -> yalmip
            DDk = sdpvar(mu,my,'full');
            eps = sdpvar(1);
            constr = (Z+U'*DDk*V+V'*DDk'*U <= eps*eye(2*n+mw+mz));
            goal = eps;
            options = sdpsettings('solver','sdpt3');
            sol = solvesdp(constr, goal, options);
            DDk = double(DDk);

        case 2  % -> 'basiclmi' from robust control toolbox
            DDk = basiclmi(Z,U,V,'Xmin,Shift');
    end
    if max(eig(Z)) <= max(eig(Z+U'*DDk*V+V'*DDk'*U))
        DDk = zeros(mu,my);
    end
end
check_LMI_DDk = max(eig(Z+U'*DDk*V+V'*DDk'*U))


% 2.3. Post-processing: auxiliary matrices
AA = A + Bu*DDk*Cy;
BBw = Bw + Bu*DDk*Dyw;
CCz = Cz + Dzu*DDk*Cy;
Dcl = Dzw + Dzu*DDk*Dyw;

Delta = -(Z+U'*DDk*V+V'*DDk'*U);
cholDelta = chol(Delta)';           % Delta = cholDelta * cholDelta'



% STEP 3. Compute BBk
% -------------------

% 3.1. Compute regular contribution
Z = [zeros(my)   , zeros(my,n), Cy         , Dyw        , zeros(my,mz);
     zeros(n,my) , -Y         , AA         , BBw        , zeros(n,mz) ;
     Cy'         , AA'        , -X         , zeros(n,mw), CCz'        ;
     Dyw'        , BBw'       , zeros(mw,n), -eye(mw)   , Dcl'        ;
     zeros(mz,my), zeros(mz,n), CCz        , Dcl        , -Gamma      ];
U = - [zeros(my,n); -eye(n); A'*X; Bw'*X; zeros(mz,n)];
BBk = Z\U;
BBk = BBk(1:my,:)';



% STEP 4. Compute CCk
% -------------------

% 4.1. Compute regular contribution
Z = [zeros(mu)   , Bu'        , zeros(mu,n), zeros(mu,mw), Dzu'       ;
     Bu          , -Y         , AA         , BBw         , zeros(n,mz);
     zeros(n,mu) , AA'        , -X         , zeros(n,mw) , CCz'       ;
     zeros(mw,mu), BBw'       , zeros(mw,n), -eye(mw)    , Dcl'       ;
     Dzu         , zeros(mz,n), CCz        , Dcl         , -Gamma     ];
U = - [zeros(mu,n); A*Y; -eye(n); zeros(mw,n); Cz*Y];
CCk = Z\U;
CCk = CCk(1:mu,:);



% STEP 5. Compute AAk
% -------------------

aux1 = cholDelta \ [-eye(n), X*A+BBk*Cy, X*Bw+BBk*Dyw, zeros(n,mz)]';
aux2 = cholDelta \ [(A*Y+Bu*CCk)', -eye(n), zeros(n,mw), (Cz*Y+Dzu*CCk)']';
AAk = - aux1'*aux2;



% STEP 6. Invert bijecive transformation of controller parameters
% ---------------------------------------------------------------

% 6.1. Compute M and N by singular value decomposition
[U,S,V] = svd(eye(n)-Y*X);
N = U*sqrt(S);
M = V*sqrt(S);


% 6.2. Reconstructing the controller
theta = [AAk,BBk;CCk,DDk] - [X*A*Y, zeros(n,my); zeros(mu,n+my)];
theta = [M, X*Bu; zeros(mu,n), eye(mu)] \ theta;
theta = theta / [N', zeros(n,my); Cy*Y, eye(my)];

Ak = theta(1:n,1:n);
Bk = theta(1:n,n+1:end);
Ck = theta(n+1:end,1:n);
Dk = theta(n+1:end,n+1:end);



% STEP 7. Verification of the result
% ----------------------------------

% 7.1. LMI in terms of the transformed controller parameters
Z = [-Y         , -eye(n)    , A*Y        , A          , Bw         , zeros(n,mz);
     -eye(n)    , -X         , zeros(n)   , X*A        , X*Bw       , zeros(n,mz);
     Y*A'       , zeros(n)   , -Y         ,-eye(n)     , zeros(n,mw), Y*Cz'      ;
     A'         , A'*X       , -eye(n)    , -X         , zeros(n,mw), Cz'        ;
     Bw'        , Bw'*X      , zeros(mw,n), zeros(mw,n), -eye(mw)   , Dzw'       ;
     zeros(mz,n), zeros(mz,n), Cz*Y       , Cz         , Dzw        , -Gamma     ];   
U = [zeros(n), eye(n)     , zeros(n,2*n) , zeros(n,mw) , zeros(n,mz);
     Bu'     , zeros(mu,n), zeros(mu,2*n), zeros(mu,mw), Dzu'       ];
V = [zeros(n,2*n) , eye(n)     , zeros(n), zeros(n,mw), zeros(n,mz) ;
     zeros(my,2*n), zeros(my,n), Cy      , Dyw        , zeros(my,mz)];
theta = [AAk,BBk;CCk,DDk]; 
check_LMI1 = max(eig(Z+U'*theta*V+V'*theta'*U));


% 7.2. LMI in terms of the actual controller parameters
Z = [-Y         , -eye(n)    , A*Y        , A          , Bw         , zeros(n,mz);
     -eye(n)    , -X         , X*A*Y      , X*A        , X*Bw       , zeros(n,mz);
     Y*A'       , Y*A'*X     , -Y         ,-eye(n)     , zeros(n,mw), Y*Cz'      ;
     A'         , A'*X       , -eye(n)    , -X         , zeros(n,mw), Cz'        ;
     Bw'        , Bw'*X      , zeros(mw,n), zeros(mw,n), -eye(mw)   , Dzw'       ;
     zeros(mz,n), zeros(mz,n), Cz*Y       , Cz         , Dzw        , -Gamma     ];   
U = [zeros(n), M'   , zeros(n,2*n) , zeros(n,mw) , zeros(n,mz);
     Bu'     , Bu'*X, zeros(mu,2*n), zeros(mu,mw), Dzu'       ];
V = [zeros(n,2*n) , N'  , zeros(n), zeros(n,mw), zeros(n,mz) ;
     zeros(my,2*n), Cy*Y, Cy      , Dyw        , zeros(my,mz)];
theta = [Ak,Bk;Ck,Dk]; 
check_LMI2 = max(eig(Z+U'*theta*V+V'*theta'*U));

check_LMI = [check_LMI1;check_LMI2]



% STEP 8. Postprocessing
% ----------------------

% 8.1. Compensation for Dyu~=0
Ak0 = Ak;   Bk0 = Bk;   Ck0 = Ck;   Dk0 = Dk;

if norm(Dyu,1) > 0,
    if norm(Dk,1) > 0,
        Myuk = eye(my)+Dyu*Dk; 
        Mkyu = eye(mu)+Dk*Dyu;
        S = svd(Myuk);
        if min(S) < sqrt(tolsing)
            error('Algebraic loop');
        else
            Ck = Mkyu\Ck;
            Dk = Mkyu\Dk;
            Ak = Ak-Bk*Dyu*Ck;
            Bk = Bk/Myuk;
        end
    else
        Ak = Ak-Bk*Dyu*Ck;
    end
end
sys_K = ss(Ak,Bk,Ck,Dk,sys_P.Ts);


% 8.2. Corresponding closed-loop system
sys_H = lft(sys_P, sys_K, mu, my);

   
