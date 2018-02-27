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

function [K,gamma]= mixedHinfsyn_MIMO8(P,W_us,my,mu,alpha,ch,options)
% mixedHinfsyn_MIMO - Controller design subject to multiple Hinf specifications on closed-loop
% transfer functions from various exogenous inputs w_j to various regulated
% outputs z_i.
%
% Inputs
%   P : generalized plant with inputs [w_j;u] and outputs [z_i;y] 
%   W_us : state space representation of the unstable output filter
%
%   my    : dimension of the measured output y
%   mu    : dimension of the control input u
%   
%   alpha : weigths determining how to treat the Hinf constraints
%           alpha_i = 0 -> ||H_i||_inf < 1
%           alpha_i > 0 -> ||H_i||_inf < gamma_i
%           objective = sum(alpha_i*gamma_i^2)
%   ch    : structure containing information concerning the channels
%           w_j->z_j
% Outputs
%   K : optimal controller
%   gamma : optimized Hinf norms: gamma_i_opt




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 6
    error('Input arguments incomplete')
else
    if(~isfield(options,'gammasolver'))
        options.gammasolver = 'lmilab';        % solver for the LMIs to compute gamma
    end
    if(~isfield(options,'controllersolver'))
        options.controllersolver = 'lmilab';   % 'basiclmi' from robust control toolbox
    end                                        % 'lmilab' (parsed with lmilab)
    if(~isfield(options,'beta'))
        options.beta = 0;
    end
    if(~isfield(options,'orderreduction'))
        options.orderreduction='off';
    end
    if(~isfield(options,'numstab'))
        options.numstab=0;
    end
end
        beta=options.beta;  %retrieve beta


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Preprocessing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%STEP 1: Retrieve plant data
%--------------------------------------------------------------------


% Step 1.1: State-space model and dimensions plant
sys_P = balreal(P);
[A,B,C,D] = ssdata(sys_P);
n = length(A(:,1));
Mz=zeros(length(ch.Out),1);
for i=1:length(Mz)
    Mz(i)=length(ch.Out(:,1));
end
Mw=zeros(length(ch.In),1);
for i=1:length(Mw)
    Mw(i)=length(ch.In(:,1));
end
mz = sum(Mz);
mw = length(B(1,:))-mu;


% Step 1.2: Subdivision according to inputs and outputs
                        Bw = B(:,1:mw);             Bu = B(:,mw+1:end);
Cz = C(1:mz,:);         Dzw = D(1:mz,1:mw);         Dzu = D(1:mz,mw+1:end);
Cy = C(mz+1:end,:);     Dyw = D(mz+1:end,1:mw);


% Step 1.3: State transformation to allow controller order reductions
% if strcmp(options.orderreduction,'on')
%     [ T,Cy_T,Dyw_T,m,index_r ] = meas_state_trans(Cy,Dyw);
%  
% else
%     T=eye(n);
      m=0;
%     Cy_T=Cy;
%     Dyw_T=Dyw;
% end
% A_T=T*A/(T);
% Bu_T=T*Bu;
% Bw_T=T*Bw;
% Cz_T=Cz/(T);
% 
% B_T=[Bw_T Bu_T];
% C_T=[Cz_T;Cy_T];
% D_T=[Dzw Dzu;Dyw_T zeros(my,mu)];
% sys_P=ss(A_T,B_T,C_T,D_T);

% Step 1.4: call controller_solver
if ~isstable(W_us)
    [K,gamma]=mixedHinfsyn_MIMO_unstab(sys_P,W_us,my,mu,m,alpha,ch,options);
else
    [K,gamma]=mixedHinfsyn_MIMO_stab(sys_P,my,mu,m,alpha,ch,options);

end

% Step 1.5: swap measurements according to inverse state transformation
if strcmp(options.orderreduction,'on')
[Ac,Bc,Cc,Dc]=ssdata(K);
for j=1:length(index_r)
    Bc(:,[j index_r(j)])=Bc(:,[index_r(j) j]);
end
K=ss(Ac,Bc,Cc,Dc);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [K,gamma]=mixedHinfsyn_MIMO_unstab(P,W_us,my,mu,m,alpha,ch,options)


beta=options.beta;
%STEP 1: Retrieve plant data
%--------------------------------------------------------------------


% Step 1.1: State-space model and dimensions plant
sys_P = P;
[A,B,C,D] = ssdata(sys_P);
n = length(A(:,1));
Mz=zeros(length(ch.Out),1);
for i=1:length(Mz)
    Mz(i)=length(ch.Out(:,1));
end
Mw=zeros(length(ch.In),1);
for i=1:length(Mw)
    Mw(i)=length(ch.In(:,1));
end
mz = sum(Mz);
mw = length(B(1,:))-mu;


% Step 1.2: Subdivision according to inputs and outputs
                        Bw = B(:,1:mw);             Bu = B(:,mw+1:end);
Cz = C(1:mz,:);         Dzw = D(1:mz,1:mw);         Dzu = D(1:mz,mw+1:end);
Cy = C(mz+1:end,:);     Dyw = D(mz+1:end,1:mw);     Dyu = D(mz+1:end,mw+1:end);

% Step 1.3: State-space model and dimensions output filter

[Ao,Bo,Co,Do]=ssdata(W_us);
no=length(Ao(:,1));

% Step 2: Build the augmented plant replicating the unstable dynamics of
% the output filter
%--------------------------------------------------------------------------

%Step 2.1: Calculate PI and gamma with yalmip
switch 'Linear' 
    case 'yalmip' %with yalmip
    %Step 2.2.1: Declare variables
    PI=sdpvar(no,n);
    Gam=sdpvar(no,my);
    
    %Step 2.2.2: Solve equation
    constr=[[PI Gam]*([A Bw;Cy Dyw])==[Bo*Cz Bo*Dzw]+[Ao*PI zeros(no,mw)],PI*Bu==Bo*Dzu];
    sol = solvesdp(constr,[],sdpsettings('solver','linprog'));
    
    
    % Step 2.2.3: Retrieve variables
    PI=double(PI);
    Gam=double(Gam);
    PI(find(isnan(PI))) = eps; %Set remaining degrees of freedom to eps
    
    case 'Linear' %solve by reshaping equalities into Ax=B
   
    %Step 2.2.1: Solve equation
    C_syl=[Bo*Cz Bo*Dzw];
    
    Lin_sol= [kron(A',eye(no))-kron(Ao,eye(n)) , kron(Cy',eye(no));
        kron(Bw',eye(no)) , kron(Dyw',eye(no))]\reshape(C_syl,numel(C_syl),1);
       
    % Step 2.2.2: Retrieve variables
    Lin_sol=reshape(Lin_sol,no,n+my);
    PI=Lin_sol(:,1:n);
    Gam=Lin_sol(:,n+1:end);
end
So=[eye(n) zeros(n,no);PI eye(no)];

% Step 2.2: Extend plant with unstable dynamics
A_hat=[A zeros(n,no);Bo*Cz Ao];
Bw_hat=[Bw;Bo*Dzw];
Bu_hat=[Bu;Bo*Dzu];
Cz_hat=[Do*Cz Co];
Dzw_hat=Do*Dzw;
Dzu_hat=Do*Dzu;
Cy_hat=[Cy zeros(my,no);-PI eye(no)];
Cy_tild=[Cy zeros(my,no);zeros(no,n) eye(no)];
Dyw_hat=[Dyw;zeros(no,mw)];



my=my+no; %The amount of measurements is extended by the order of the unstable weight
z=PI*Bu-Bo*Dzu;%For later use


% Step 3: Pre-processing of matrices to account for multiple objectives
%--------------------------------------------------------------------------


% Step 3.1: Regulated outputs
alpha = alpha(:);
noc=length(ch.In);
Iz2 = cumsum(Mz);
Iz1 = [1;Iz2(1:end-1)+1];

% Step 3.2: External inputs
Iw2 = cumsum(Mw);
Iw1 = [1;Iw2(1:end-1)+1];
mw=noc;

%Step 3.3: Duplicate matrices

Bw_hat_dub=Bw_hat;
Dzw_hat_dub=Dzw_hat;
Dyw_hat_dub=Dyw_hat;
Bw_dub=Bw;

Bw=[];
Bw_hat=[];
Dzw_hat=[];
Dyw_hat=[];

for i=1:noc
j=find(ch.In{i});
Bw_hat  = [Bw_hat Bw_hat_dub(:,j)];
Bw  = [Bw Bw_dub(:,j)];
Dzw_hat  = [Dzw_hat Dzw_hat_dub(:,j)];
Dyw_hat  = [Dyw_hat Dyw_hat_dub(:,j)];
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Controller design
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% STEP 0: Settings
% ----------------

if strcmp(options.gammasolver, 'lmilab')
    LMIparser = 'lmilab';                      % 2: LMI Lab
else
    LMIparser = 'yalmip';                      % 1: yalmip
end



% STEP 1: SDP to compute the optimal gamma_i
% ------------------------------------------

% Step 1.1: Preprocessing: null-spaces and auxiliary matrices
V = null([Cy_tild, Dyw_hat]);
id_V=[1:n-m n+no+1:length(V(:,1))];%Remove null-rows
V=V(id_V,:);
V = blkdiag(V, eye(mz));

W = null([Bu_hat', Dzu_hat']);
W = blkdiag(W, eye(mw));

Ab=A_hat+beta*eye(n+no);
Sigma=Do*Cz+Co*PI;

switch LMIparser
  
    case 'yalmip' %parse with YALMIP, solve with Mosek

        % Step 1b.1: Declaring the variables
        X11 = sdpvar(n-m,n-m,'symmetric');
        X21= sdpvar(m, n-m);
        X31=sdpvar(no,n-m);
        X= [X11; X21];
        Y = sdpvar(n+no,n+no,'symmetric');
        gamma = sdpvar(noc,1,'full');
        Gamma_w = sdpvar(mw);
        Gamma_z = sdpvar(mz);
        DDzw = sdpvar(mz,mw,'full');
        
        for i = 1:noc
            if alpha(i) == 0
                gamma(i) = 1;
            end
            
            Gamma_z(Iz1(i):Iz2(i),Iz1(i):Iz2(i)) = gamma(i) * eye(Mz(i)); % eye(Mz(i))
            Gamma_w(Iw1(i):Iw2(i),Iw1(i):Iw2(i)) = gamma(i) * eye(Mw(i)); % eye(Mw(i))
            DDzw(Iz1(i):Iz2(i),Iw1(i):Iw2(i)) = Dzw_hat(Iz1(i):Iz2(i),Iw1(i):Iw2(i));
        end
        % Step 1b.2: Specifying the LMIs
        Zx = [A(:,1:n-m)'*X+X'*A(:,1:n-m), X'*Bw    , Sigma(:,1:n-m)'   ;
            Bw'*X   , -Gamma_w, DDzw'  ;
            Sigma(:,1:n-m)      , DDzw     , -Gamma_z];
        
        Zy = [Ab*So*Y*So'+So*Y*So'*Ab', So*Y*So'*Cz_hat' , Bw_hat      ;
            Cz_hat*So*Y*So'    , -Gamma_z, DDzw     ;
            Bw_hat'     , DDzw'  , -Gamma_w];
        
        I_XY_pos=[eye(n-m) zeros(n-m,m+no)];
        XY_pos=[Y I_XY_pos';I_XY_pos X11];
        constr = [XY_pos >= 10e-8 , V'*Zx*V <= eps*eye(length(V(1,:))) , W'*Zy*W <= eps*eye(length(W(1,:)))];
        
        
        % Step 1.4: Specifying the objective
        goal = trace(alpha'*gamma);
        
        
        % Step 1.5: Solving the SDP
        sol = optimize(constr, goal, sdpsettings('solver',options.gammasolver))
        check_LMIs = {'P pos' , min(eig(double(XY_pos)));
            'Zx neg', max(eig(double(V'*Zx*V)))            ;
            'Zy neg', max(eig(double(W'*Zy*W)))            }
        
        % Step 1.6: Retrieve variables
        X11 = value(X11);
        X21 = [double(X21);double(X31)];
        X21(find(isnan(X21))) = eps;
        Y = value(Y);
        XY_pos=double(XY_pos);
        gamma =  (value(gamma)); % value(gamma); % sqrt(value(gamma));
        Gamma_z = value(Gamma_z);
        Gamma_w = value(Gamma_w);
        DDzw = value(DDzw);
        
    case 'lmilab' %Parse with LMILab, solve with LMILab

        setlmis([])
        
        % Step 1.2: Declaring the variables
        [X11,~,sX11] = lmivar(1,[n-m,1]);
        [X21,~,sX21] = lmivar(2,[m,n-m]);
        [X31,~,sX31] = lmivar(2,[no,n-m]);
        [X,~,sX]=lmivar(3,[sX11;sX21]);
        [Y,mvar,~] = lmivar(1,[n+no,1]);
        struc_gam=[mvar+1:mvar+1+noc];
        Gamma0_z = zeros(mz);
        Gamma0_w = zeros(mw);
        strucGamma_w = zeros(mw);
        strucGamma_z = zeros(mz);
        DDzw_fix = zeros(mz,mw);
        strucDDzw_var=zeros(mz,mw);
        for i = 1:noc
            if alpha(i) == 0
                Gamma0_z(Iz1(i):Iz2(i),Iz1(i):Iz2(i)) = eye(Mz(i));
                Gamma0_w(Iw1(i):Iw2(i),Iw1(i):Iw2(i)) = eye(Mw(i));
            else
                strucGamma_z(Iz1(i):Iz2(i),Iz1(i):Iz2(i)) = ((mvar+1)/2) *eye(Mz(i));
                strucGamma_w(Iw1(i):Iw2(i),Iw1(i):Iw2(i)) = ((mvar+1)/2) *eye(Mw(i));
                mvar = mvar+1;
            end
            temp = reshape(mvar+(1:Mz(i)*sum(Mz(i+1:end)))', Mz(i), sum(Mz(i+1:end)));
            strucGamma_z(Iz1(i):Iz2(i),Iz2(i)+1:end) = temp;
            mvar = mvar + Mz(i)*sum(Mz(i+1:end));
            
            temp = reshape(mvar+(1:Mw(i)*sum(Mw(i+1:end)))', Mw(i), sum(Mw(i+1:end)));
            strucGamma_w(Iw1(i):Iw2(i),Iw2(i)+1:end) = temp;
            mvar = mvar + Mw(i)*sum(Mw(i+1:end));
            
            DDzw_fix(Iz1(i):Iz2(i),Iw1(i):Iw2(i)) = Dzw_hat(Iz1(i):Iz2(i),Iw1(i):Iw2(i));
            strucDDzw_var(Iz1(i):Iz2(i),Iw1(i):Iw2(i))=1;
        end
        
        [strucGamma_z]=strucGamma_z+strucGamma_z';
        [Gamma_z,~,sGz] = lmivar(3, strucGamma_z);
        
        [strucGamma_w]=strucGamma_w+strucGamma_w';
        [Gamma_w,~,sGw] = lmivar(3, strucGamma_w);
        
        strucDDzw_comp=find(strucDDzw_var);
        strucDDzw_var=zeros(mz,mw);
        strucDDzw_diff=setdiff([1:numel(strucDDzw_var)],strucDDzw_comp);
        strucDDzw_var(strucDDzw_diff)=[mvar+(1:numel(strucDDzw_diff))];
        
        [DDzw_var,~,sDDzwv]=lmivar(3,strucDDzw_var);
        
        % Step 1.3: Specifying the LMIs
      
        IO=[eye(n-m) zeros(n-m,mw)];
        OI=[zeros(mw,n-m) eye(mw)];
        
        Zx_neg=newlmi;
        lmiterm([Zx_neg,1,1,X],[A(:,1:n-m) Bw]',IO,'s');
        lmiterm([Zx_neg,1,1,Gamma_w],-OI',OI);
        lmiterm([Zx_neg,1,1,0],-OI'*Gamma0_w*OI);
        lmiterm([Zx_neg,2,1,0],[Sigma(:,1:n-m) DDzw_fix]);
        lmiterm([Zx_neg,2,1,DDzw_var],1,[zeros(n-m,mw); eye(mw)]');
        lmiterm([Zx_neg,2,2,Gamma_z],-1,1);
        lmiterm([Zx_neg,2,2,0],-Gamma0_z);
        lmiterm([Zx_neg,0,0,0],V);
        
        IO = [eye(n+no), zeros(n+no,mz)];
        OI = [zeros(mz,n+no), eye(mz)];
        
        Zy_neg=newlmi;
        lmiterm([Zy_neg,1,1,Y],[Ab*So;Cz_hat*So],So'*IO,'s');
        lmiterm([Zy_neg,1,1,Gamma_z],-OI',OI);
        lmiterm([Zy_neg,1,1,0],-OI'*Gamma0_z*OI);
        lmiterm([Zy_neg,1,2,0],[Bw_hat; DDzw_fix]);
        lmiterm([Zy_neg,1,2,DDzw_var],[zeros(n+no,mw); eye(mw)],1);
        lmiterm([Zy_neg,2,2,Gamma_w],-1,1);
        lmiterm([Zy_neg,2,2,0],-Gamma0_w);
        lmiterm([Zy_neg,0,0,0],W);
        
      
        Zxy_pos=newlmi;
        I_XY_pos=[eye(n-m) zeros(n-m,m+no)];
        lmiterm([-Zxy_pos 1 1 Y],1,1);
        lmiterm([-Zxy_pos 1 2 0],I_XY_pos');
        lmiterm([-Zxy_pos 2 2 X11],1,1);
        lmiterm([Zxy_pos  1 1 0],1e-3);
        
        LMIs=getlmis;
        
        % Step 1.4: Specifying the objective
        c = zeros(decnbr(LMIs),1);
        diagGamma_z = diag(decinfo(LMIs,Gamma_z));
        diagGamma_z = diagGamma_z(Iz1);
        c(diagGamma_z(diagGamma_z~=0)) = alpha((diagGamma_z~=0));
        
        % Step 1.5: Solving the SDP
        
        [fopt, xopt] = mincx(LMIs, c);
        
        LMIsopt = evallmi(LMIs, xopt);
        [Zx_neg, trash] = showlmi(LMIsopt, Zx_neg);
        [Zy_neg, trash] = showlmi(LMIsopt, Zy_neg);
        [trash, Zxy_pos] = showlmi(LMIsopt, Zxy_pos);
        check_LMIs = {'XY pos' , min(eig(Zxy_pos)) ;
            'X neg', max(eig(Zx_neg));
            'Y neg', max(eig(Zy_neg))}
        % Step 1.6: Retrieve variables
        X11 = dec2mat(LMIs, xopt,X11);
        X21 = dec2mat(LMIs, xopt,X21);
        X31 = dec2mat(LMIs, xopt, X31);
        X21 = [double(X21);double(X31)];
        X21(find(isnan(X21))) = eps;
        X = dec2mat(LMIs, xopt,X);
        Y = dec2mat(LMIs, xopt,Y);
        Gamma_z = dec2mat(LMIs, xopt,Gamma_z)+Gamma0_z;
        Gamma_w = dec2mat(LMIs, xopt,Gamma_w)+Gamma0_w;
        gamma=(diag(Gamma_z));
        XY_pos=[Y I_XY_pos';I_XY_pos X11];
        DDzw=dec2mat(LMIs,xopt,DDzw_var)+DDzw_fix;

end

% Step 1bis: Improve numerical conditioning
% -------------------------------------------
if options.numstab %Force eigenvalues of P further away from the imaginary axis
    gamma=gamma*1.01;
switch LMIparser
    
    case 'yalmip' %Parse with YALMIP, solve with SeDuMi
        
        % Step 1b.1: Declaring the variables
        X11 = sdpvar(n-m,n-m,'symmetric');
        X21= sdpvar(m, n-m);
        X31=sdpvar(no,n-m);
        X= [X11; X21];
        Y = sdpvar(n+no,n+no,'symmetric');
        Gamma_w = sdpvar(mw);
        Gamma_z = sdpvar(mz);
        DDzw = sdpvar(mz,mw,'full');
        t=sdpvar(1);
        for i = 1:noc
            if alpha(i) == 0
                gamma(i) = 1;
            end
            
            Gamma_z(Iz1(i):Iz2(i),Iz1(i):Iz2(i)) = gamma(i) * eye(Mz(i)); % eye(Mz(i))
            Gamma_w(Iw1(i):Iw2(i),Iw1(i):Iw2(i)) = gamma(i) * eye(Mw(i)); % eye(Mw(i))
            DDzw(Iz1(i):Iz2(i),Iw1(i):Iw2(i)) = Dzw_hat(Iz1(i):Iz2(i),Iw1(i):Iw2(i));
        end
        
        % Step 1b.2: Specifying the LMIs
       
        Zx = [A(:,1:n-m)'*X+X'*A(:,1:n-m), X'*Bw    , Sigma(:,1:n-m)'   ;
            Bw'*X   , -Gamma_w, DDzw'  ;
            Sigma(:,1:n-m)      , DDzw     , -Gamma_z];
        
        Zy = [Ab*So*Y*So'+So*Y*So'*Ab', So*Y*So'*Cz_hat' , Bw_hat      ;
            Cz_hat*So*Y*So'    , -Gamma_z, DDzw     ;
            Bw_hat'     , DDzw'  , -Gamma_w];
        
        I_XY_pos=[eye(n-m) zeros(n-m,m+no)];
        XY_pos=[Y I_XY_pos';I_XY_pos X11];
        constr = [XY_pos >= t , V'*Zx*V <= -1e-8 , W'*Zy*W <= -1e-8];
        
        
        % Step 1b.3: Specifying the objective
        goal = -t;
        
        
        % Step 1b.4: Solving the SDP
        sol = optimize(constr, goal, sdpsettings('solver',options.gammasolver))
        check_LMIs = {'P pos' , min(eig(double(XY_pos)));
            'Zx neg', max(eig(double(V'*Zx*V)))            ;
            'Zy neg', max(eig(double(W'*Zy*W)))            ;
            't', double(t)}
        
        % Step 1b.5: Retrieve variables
        X11 = double(X11);
        X21 = [double(X21);double(X31)];
        X21(find(isnan(X21))) = eps;
        Y = double(Y);
        XY_pos=double(XY_pos);
        Gamma_z = double(Gamma_z);
        Gamma_w = double(Gamma_w);
        DDzw = double(DDzw);
        t=double(t);
        
    case 'lmilab' %Parse with LMILab, solve with LMILab
        
         setlmis([])
        
        % Step 1.2b : Declaring the variables
        [X11,~,sX11] = lmivar(1,[n-m,1]);
        [X21,~,sX21] = lmivar(2,[m,n-m]);
        [X31,~,sX31] = lmivar(2,[no,n-m]);
        [X,~,sX]=lmivar(3,[sX11;sX21]);
        [t,mt]=lmivar(1,[1,1]);
        [Y,mvar] = lmivar(1,[n+no,1]);
        Gamma0_z = zeros(mz);
        Gamma0_w = zeros(mw);
        strucGamma_w = zeros(mw);
        strucGamma_z = zeros(mz);
        DDzw_fix = zeros(mz,mw);
        strucDDzw_var=zeros(mz,mw);
        for i = 1:noc
            if alpha(i) == 0
                Gamma0_z(Iz1(i):Iz2(i),Iz1(i):Iz2(i)) = eye(Mz(i));
                Gamma0_w(Iw1(i):Iw2(i),Iw1(i):Iw2(i)) = eye(Mw(i));
            else
                Gamma0_z(Iz1(i):Iz2(i),Iz1(i):Iz2(i)) = gamma(i);
                Gamma0_w(Iw1(i):Iw2(i),Iw1(i):Iw2(i)) = gamma(i);
            end
            temp = reshape(mvar+(1:Mz(i)*sum(Mz(i+1:end)))', Mz(i), sum(Mz(i+1:end)));
            strucGamma_z(Iz1(i):Iz2(i),Iz2(i)+1:end) = temp;
            mvar = mvar + Mz(i)*sum(Mz(i+1:end));
            
            temp = reshape(mvar+(1:Mw(i)*sum(Mw(i+1:end)))', Mw(i), sum(Mw(i+1:end)));
            strucGamma_w(Iw1(i):Iw2(i),Iw2(i)+1:end) = temp;
            mvar = mvar + Mw(i)*sum(Mw(i+1:end));
            
            DDzw_fix(Iz1(i):Iz2(i),Iw1(i):Iw2(i)) = Dzw_hat(Iz1(i):Iz2(i),Iw1(i):Iw2(i));
            strucDDzw_var(Iz1(i):Iz2(i),Iw1(i):Iw2(i))=1;
        end
        
        [strucGamma_z]=strucGamma_z+strucGamma_z';
        [Gamma_z,~,sGz] = lmivar(3, strucGamma_z);
        
        [strucGamma_w]=strucGamma_w+strucGamma_w';
        [Gamma_w,~,sGw] = lmivar(3, strucGamma_w);
        
        strucDDzw_comp=find(strucDDzw_var);
        strucDDzw_var=zeros(mz,mw);
        strucDDzw_diff=setdiff([1:numel(strucDDzw_var)],strucDDzw_comp);
        strucDDzw_var(strucDDzw_diff)=[mvar+(1:numel(strucDDzw_diff))];
        
        if sum(strucDDzw_var(:))~=0
            [DDzw_var,~,sDDzwv]=lmivar(3,strucDDzw_var);
        end        
        
        % Step 1.3b : Specifying the LMIs
        
        IO=[eye(n-m) zeros(n-m,mw)];
        OI=[zeros(mw,n-m) eye(mw)];
        
        Zx_neg=newlmi;
        lmiterm([Zx_neg,1,1,X],[A(:,1:n-m) Bw]',IO,'s');
        lmiterm([Zx_neg,1,1,Gamma_w],-OI',OI);
        lmiterm([Zx_neg,1,1,0],-OI'*Gamma0_w*OI);
        lmiterm([Zx_neg,2,1,0],[Sigma(:,1:n-m) DDzw_fix]);
        if sum(strucDDzw_var(:))~=0; lmiterm([Zx_neg,2,1,DDzw_var],1,[zeros(n-m,mw); eye(mw)]'); end;
        lmiterm([Zx_neg,2,2,Gamma_z],-1,1);
        lmiterm([Zx_neg,2,2,0],-Gamma0_z);
        lmiterm([Zx_neg,0,0,0],V);
        
        IO = [eye(n+no), zeros(n+no,mz)];
        OI = [zeros(mz,n+no), eye(mz)];
        
        Zy_neg=newlmi;
        lmiterm([Zy_neg,1,1,Y],[Ab*So;Cz_hat*So],So'*IO,'s');
        lmiterm([Zy_neg,1,1,Gamma_z],-OI',OI);
        lmiterm([Zy_neg,1,1,0],-OI'*Gamma0_z*OI);
        lmiterm([Zy_neg,1,2,0],[Bw_hat; DDzw_fix]);
        if sum(strucDDzw_var(:))~=0; lmiterm([Zy_neg,1,2,DDzw_var],[zeros(n+no,mw); eye(mw)],1); end;
        lmiterm([Zy_neg,2,2,Gamma_w],-1,1);
        lmiterm([Zy_neg,2,2,0],-Gamma0_w);
        lmiterm([Zy_neg,0,0,0],W);
        
        
        Zxy_pos=newlmi;
        I_XY_pos=[eye(n-m) zeros(n-m,m+no)];
        lmiterm([-Zxy_pos 1 1 Y],1,1);
        lmiterm([-Zxy_pos 1 2 0],I_XY_pos');
        lmiterm([-Zxy_pos 2 2 X11],1,1);
        lmiterm([Zxy_pos 1 1 t],1,1);
        
        t_pos=newlmi;
        lmiterm([-t_pos 1 1 t],1,1);
        lmiterm([-t_pos 2 2 t],1,1);

        
        LMIs=getlmis;
        
        % Step 1.4b : Specifying the objective
        c = zeros(decnbr(LMIs),1);
        c(mt) = -1;
        
        % Step 1.5b : Solving the SDP
        [fopt, xopt] = mincx(LMIs, c);
        
        LMIsopt = evallmi(LMIs, xopt);
        [Zx_neg, trash] = showlmi(LMIsopt, Zx_neg);
        [Zy_neg, trash] = showlmi(LMIsopt, Zy_neg);
        [trash, Zxy_pos] = showlmi(LMIsopt, Zxy_pos);
%        [trash, t_pos] = showlmi(LMIsopt, t_pos);

        check_LMIs = {'XY pos' , min(eig(Zxy_pos)) ;
            'X neg', max(eig(Zx_neg));
            'Y neg', max(eig(Zy_neg));
            't_pos', min(eig(t_pos))
 }
        
        % Step 1.6b : Retrieve variables
        X11 = dec2mat(LMIs, xopt,X11);
        X21 = dec2mat(LMIs, xopt,X21);
        X31 = dec2mat(LMIs, xopt, X31);
        X21 = [double(X21);double(X31)];
        X21(find(isnan(X21))) = eps;
        Y = dec2mat(LMIs, xopt,Y);
        Gamma_z = dec2mat(LMIs, xopt,Gamma_z)+Gamma0_z;
        Gamma_w = dec2mat(LMIs, xopt,Gamma_w)+Gamma0_w;
        XY_pos=[Y I_XY_pos';I_XY_pos X11];
        DDzw=dec2mat(LMIs,xopt,DDzw_var)+DDzw_fix;
end
end
% STEP 2: Reconstructing the controller
% -------------------------------------

%Step 2.1: Calculate M,N
[U,S,V] = svd(eye(n-m)-[X11 X21']*[Y(1:n-m,1:n-m);Y(n-m+1:n+no,1:n-m)]);
N1 = V*sqrt(S);
M = U*sqrt(S);
N2=-[Y(n-m+1:n+no,1:n-m) Y(n-m+1:end,n-m+1:end)]*[X11;X21]/(M');

%Step 3.2: Calculate Lyapunov matrix P trough reverse congruence transform
PI_c=[Y [eye(n-m);zeros(m+no,n-m)]; N1' N2' zeros(n-m)];
P=PI_c'\XY_pos/PI_c;
P=blkdiag(So,eye(n-m))'\P/blkdiag(So,eye(n-m));

% Step 2.2: Retrieving the controller
% switch 1
%     case 1 %Controller parameters
nc = n-m;
Aa = blkdiag(A_hat+beta*eye(n+no),beta*eye(nc));
Bwa = [Bw_hat; zeros(nc,mw)];
Bua = [zeros(n+no,nc), Bu_hat; eye(nc), zeros(nc,mu)];
Cza = [Cz_hat, zeros(mz,nc)];
Dzua = [zeros(mz,nc), Dzu_hat];

%reverse measurements, replace by better code
Cy_hat=[-PI eye(no);Cy zeros(my-no,no)];
Dyw_hat=[zeros(no,size(ch.In{1}',2));Dyw];
Dyw_hat_dub=Dyw_hat;
Dyw_hat=[];
for i=1:noc
    j=find(ch.In{i});
    Dyw_hat  = [Dyw_hat Dyw_hat_dub(:,j)];
end

Cya = [zeros(nc,n+no), eye(nc); Cy_hat, zeros(my,nc)];
Dywa = [zeros(nc,mw); Dyw_hat];
Z0 = [Aa'*P+P*Aa, P*Bwa     , Cza';
    Bwa'*P      , -Gamma_w, DDzw'  ;
    Cza     , DDzw     , -Gamma_z];
L = [P*Bua; zeros(mw,nc+mu); Dzua]';
R = [Cya, Dywa, zeros(nc+my,mz)];

%     case 2  % -> transformed controller parameters. BUGGED, use exact controller parameters (case 1) instead
%         nc = n-m;
%         
%         X1=[X11 X21']
%         PI_I=[eye(n-m);zeros(m+no,n-m)]
%         
%         Z0 = [A_hat*Y+Y*A_hat', A_hat*PI_I       , Bw_hat      , Y*Cz_hat' ;
%             PI_I'*A_hat'      , X1*A_hat*PI_I+PI_I'*A_hat'*X1', X1*Bw_hat    , PI_I'*Cz_hat'   ;
%             Bw_hat'     , Bw_hat'*X1'   , -Gamma_w, DDzw'  ;
%             Cz_hat*Y    , Cz_hat*PI_I      , DDzw     , -Gamma_z];
%         L = [zeros(nc,n+no), eye(nc)     , zeros(nc,mw) , zeros(nc,mz);
%             Bu_hat'     , zeros(mu,nc), zeros(mu,mw), Dzu_hat'       ];
%         R = [eye(n+no)     , zeros(n+no,nc), zeros(n+no,mw), zeros(n+no,mz) ;
%             zeros(my,n+no), Cy_hat*PI_I      , Dyw_hat        , zeros(my,mz)];
% end

switch options.controllersolver
    case 'basiclmi' %Basiclmi
        theta = basiclmi(Z0,L,R,'Xmin');
    case 'lmilab' %LMIlab
        setlmis([]);
        theta=lmivar(2,[nc+mu nc+my]);
        theta_calc=newlmi;
        lmiterm([theta_calc 1 1 0],Z0);
        lmiterm([theta_calc 1 1 theta],L',R,'s')
        LMIs=getlmis;
        [fopt, xopt] = feasp(LMIs);
    
        theta = dec2mat(LMIs, xopt, theta);
        
    otherwise       % -> yalmip + SDPsolver
        theta = sdpvar(nc+mu,nc+my,'full');
        t = sdpvar(1);
        constr = [Z0+L'*theta*R+R'*theta'*L <= t*eye(n+no+nc+mw+mz)];
        options = sdpsettings('solver',options.controllersolver);
        goal = t;
        sol = solvesdp(constr, goal, options);
        theta = double(theta);
        check_t = value(t)
end
check_LMI = max(eig(Z0+L'*theta*R+R'*theta'*L))


% Step 3: Recombining the controller
%--------------------------------------

%Step 3.1: partition theta
Ac_22 = theta(1:nc,1:nc);
Ac_21 = theta(1:nc,nc+1:nc+no);
Bc_2= theta(1:nc,nc+no+1:end);
Cc_2 = theta(nc+1:end,1:nc);
Cc_1 =theta(nc+1:end,nc+1:nc+no);
Dc0=theta(nc+1:end,nc+no+1:end);

% Ac_hat=Ac_22;
% Bc_hat=[Ac_21 Bc_2];
% Cc_hat=[Cc_2];
% Dc_hat=[Cc_1 Dc];
% 
% %step 3.2: Calculate predetermined controller parts
% Acl_test=[A_hat+Bu_hat*[Cc_1 Dc]*Cy_hat Bu_hat*[Cc_2];[Ac_21 Bc_2]*Cy_hat Ac_22];
% if max(real(eig(Acl_test)))>=0
%     warning('augmented closed loop has not been stabilized. Controller design failed.')
% end
Ac_11=Ao-z*Cc_1;
Ac_12=-z*Cc_2;
Bc_1=(Gam-z*Dc0);

%step 3.3: Combine
Ac0=[Ac_11 Ac_12;Ac_21 Ac_22];
Bc0=[Bc_1;Bc_2];
Cc0=[Cc_1 Cc_2];

if norm(Dyu,1) > 0,
        if norm(Dk0,1) > 0,
            Myuk = eye(dims.my)+Dyu*Dc0; 
            Mkyu = eye(dims.mu)+Dc0*Dyu;
            if svds(Myuk,1,'smallest') < 1e-6
                error('The controller could not be reconstructed due to a singularity. This is usually caused by an algebraic loop in the generalized plant.');
            else
                Ck = Mkyu\Cc0;
                Dk = Mkyu\Dc0;
                Ak = Ak0-Bk0*Dyu*Ck;
                Bk = Bk0/Myuk;
            end
        else
            Ak = Ac0-Bc0*Dyu*Cc0;
        end
    else
         Ak = Ac0; Bk = Bc0; Ck = Cc0; Dk = Dc0;
end

K=ss(Ak,Bk,Ck,Dk);
if ~isstable(lft(sys_P,K));
    warning('closed loop has not been stabilized. Controller design failed.')
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [K,gamma]=mixedHinfsyn_MIMO_stab(P,my,mu,m,alpha,ch,options)
beta=options.beta;
%STEP 1: Retrieve plant data
%--------------------------------------------------------------------


% Step 1.1: State-space model and dimensions plant
sys_P = P;
[A,B,C,D] = ssdata(sys_P);
n = length(A(:,1));
Mz=zeros(length(ch.Out),1);
for i=1:length(Mz)
    Mz(i)=length(ch.Out(:,1));
end
Mw=zeros(length(ch.In),1);
for i=1:length(Mw)
    Mw(i)=length(ch.In(:,1));
end
mz = sum(Mz);
mw = length(B(1,:))-mu;


% Step 1.2: Subdivision according to inputs and outputs
                        Bw = B(:,1:mw);             Bu = B(:,mw+1:end);
Cz = C(1:mz,:);         Dzw = D(1:mz,1:mw);         Dzu = D(1:mz,mw+1:end);
Cy = C(mz+1:end,:);     Dyw = D(mz+1:end,1:mw);


% Step 2: Pre-processing of matrices to account for multiple objectives
%--------------------------------------------------------------------------


% Step 2.1: Regulated outputs
alpha = alpha(:);
noc=length(ch.In);
Iz2 = cumsum(Mz);
Iz1 = [1;Iz2(1:end-1)+1];

% Step 2.2: External inputs
Iw2 = cumsum(Mw);
Iw1 = [1;Iw2(1:end-1)+1];
mw=noc;

%Step 2.3: Duplicate matrices

Bw_dub=Bw;
Dzw_dub=Dzw;
Dyw_dub=Dyw;

Bw=[];
Dzw=[];
Dyw=[];

for i=1:noc
j=find(ch.In{i});
Bw       = [Bw Bw_dub(:,j)];
Dzw  = [Dzw Dzw_dub(:,j)];
Dyw  = [Dyw Dyw_dub(:,j)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Controller design
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% STEP 0: Settings
% ----------------

if strcmp(options.gammasolver, 'lmilab')
    LMIparser = 'lmilab';                      % 2: LMI Lab
else
    LMIparser = 'yalmip';                      % 1: yalmip
end



% STEP 1: SDP to compute the optimal gamma_i
% ------------------------------------------

% Step 1.1: Preprocessing: null-spaces and auxiliary matrices
V = null([Cy, Dyw]);
id_V=[1:n-m n+1:length(V(:,1))];%Remove null-rows
V=V(id_V,:);
V = blkdiag(V, eye(mz));

W = null([Bu', Dzu']);
W = blkdiag(W, eye(mw));

Ab=A+beta*eye(n);

switch LMIparser
  
    case 'yalmip' %parse with YALMIP, solve with Mosek

        % Step 1b.1: Declaring the variables
        X11 = sdpvar(n-m,n-m,'symmetric');
        X21= sdpvar(m, n-m);
        X= [X11; X21];
        Y = sdpvar(n,n,'symmetric');
        gamma = sdpvar(noc,1,'full');
        Gamma_w = sdpvar(mw);
        Gamma_z = sdpvar(mz);
        DDzw = sdpvar(mz,mw,'full');
        
        for i = 1:noc
            if alpha(i) == 0
                gamma(i) = 1;
            end
            
            Gamma_z(Iz1(i):Iz2(i),Iz1(i):Iz2(i)) = gamma(i) * eye(Mz(i)); % eye(Mz(i))
            Gamma_w(Iw1(i):Iw2(i),Iw1(i):Iw2(i)) = gamma(i) * eye(Mw(i)); % eye(Mw(i))
            DDzw(Iz1(i):Iz2(i),Iw1(i):Iw2(i)) = Dzw(Iz1(i):Iz2(i),Iw1(i):Iw2(i));
        end
        % Step 1b.2: Specifying the LMIs
        Zx = [Ab(:,1:n-m)'*X+X'*Ab(:,1:n-m), X'*Bw    , Cz(:,1:n-m)'   ;
            Bw'*X   , -Gamma_w, DDzw'  ;
            Cz(:,1:n-m)      , DDzw     , -Gamma_z];
        
        Zy = [Ab*Y'+Y*Ab', Y*Cz' , Bw      ;
            Cz*Y    , -Gamma_z, DDzw     ;
            Bw'     , DDzw'  , -Gamma_w];
        
        I_XY_pos=[eye(n-m) zeros(n-m,m)];
        XY_pos=[Y I_XY_pos';I_XY_pos X11];
        constr = [XY_pos >= 10e-8 , V'*Zx*V <= eps*eye(length(V(1,:))) , W'*Zy*W <= eps*eye(length(W(1,:)))];
        
        
        % Step 1.4: Specifying the objective
        goal = trace(alpha'*gamma);
        
        
        % Step 1.5: Solving the SDP
        sol = optimize(constr, goal, sdpsettings('solver',options.gammasolver))
        check_LMIs = {'P pos' , min(eig(double(XY_pos)));
            'Zx neg', max(eig(double(V'*Zx*V)))            ;
            'Zy neg', max(eig(double(W'*Zy*W)))            }
        
        % Step 1.6: Retrieve variables
        X11 = value(X11);
        X21 = [double(X21)];
        X21(find(isnan(X21))) = eps;
        Y = value(Y);
        XY_pos=double(XY_pos);
        gamma =  (value(gamma)); % value(gamma); % sqrt(value(gamma));
        Gamma_z = value(Gamma_z);
        Gamma_w = value(Gamma_w);
        DDzw = value(DDzw);
        
    case 'lmilab' %Parse with LMILab, solve with LMILab
        
        setlmis([])
        
        % Step 1.2: Declaring the variables
        [X11,~,sX11] = lmivar(1,[n-m,1]);
        [X21,~,sX21] = lmivar(2,[m,n-m]);
        [X,~,sX]=lmivar(3,[sX11;sX21]);
        [Y,mvar,sY] = lmivar(1,[n,1]);
        Gamma0_z = zeros(mz);
        Gamma0_w = zeros(mw);
        strucGamma_w = zeros(mw);
        strucGamma_z = zeros(mz);
        DDzw_fix = zeros(mz,mw);
        strucDDzw_var=zeros(mz,mw);
        for i = 1:noc
            if alpha(i) == 0
                Gamma0_z(Iz1(i):Iz2(i),Iz1(i):Iz2(i)) = eye(Mz(i));
                Gamma0_w(Iw1(i):Iw2(i),Iw1(i):Iw2(i)) = eye(Mw(i));
            else
                strucGamma_z(Iz1(i):Iz2(i),Iz1(i):Iz2(i)) = ((mvar+1)/2) *eye(Mz(i));
                strucGamma_w(Iw1(i):Iw2(i),Iw1(i):Iw2(i)) = ((mvar+1)/2) *eye(Mw(i));
                mvar = mvar+1;
            end
            temp = reshape(mvar+(1:Mz(i)*sum(Mz(i+1:end)))', Mz(i), sum(Mz(i+1:end)));
            strucGamma_z(Iz1(i):Iz2(i),Iz2(i)+1:end) = temp;
            mvar = mvar + Mz(i)*sum(Mz(i+1:end));
            
            temp = reshape(mvar+(1:Mw(i)*sum(Mw(i+1:end)))', Mw(i), sum(Mw(i+1:end)));
            strucGamma_w(Iw1(i):Iw2(i),Iw2(i)+1:end) = temp;
            mvar = mvar + Mw(i)*sum(Mw(i+1:end));
            
            DDzw_fix(Iz1(i):Iz2(i),Iw1(i):Iw2(i)) = Dzw(Iz1(i):Iz2(i),Iw1(i):Iw2(i));
            strucDDzw_var(Iz1(i):Iz2(i),Iw1(i):Iw2(i))=1;
        end
        
        [strucGamma_z]=strucGamma_z+strucGamma_z';
        [Gamma_z,~,sGz] = lmivar(3, strucGamma_z);
        
        [strucGamma_w]=strucGamma_w+strucGamma_w';
        [Gamma_w,~,sGw] = lmivar(3, strucGamma_w);
        
        strucDDzw_comp=find(strucDDzw_var);
        strucDDzw_var=zeros(mz,mw);
        strucDDzw_diff=setdiff([1:numel(strucDDzw_var)],strucDDzw_comp);
        strucDDzw_var(strucDDzw_diff)=[mvar+(1:numel(strucDDzw_diff))];
        
        if sum(strucDDzw_var(:))~=0
            [DDzw_var,~,sDDzwv]=lmivar(3,strucDDzw_var);
        end
        
        % Step 1.3: Specifying the LMIs
      
        IO=[eye(n-m) zeros(n-m,mw)];
        OI=[zeros(mw,n-m) eye(mw)];
        
        Zx_neg=newlmi;
        lmiterm([Zx_neg,1,1,X],[Ab(:,1:n-m) Bw]',IO,'s');
        lmiterm([Zx_neg,1,1,Gamma_w],-OI',OI);
        lmiterm([Zx_neg,1,1,0],-OI'*Gamma0_w*OI);
        lmiterm([Zx_neg,2,1,0],[Cz(:,1:n-m) DDzw_fix]);
        if sum(strucDDzw_var(:))~=0; lmiterm([Zx_neg,2,1,DDzw_var],1,[zeros(n-m,mw); eye(mw)]'); end;
        lmiterm([Zx_neg,2,2,Gamma_z],-1,1);
        lmiterm([Zx_neg,2,2,0],-Gamma0_z);
        lmiterm([Zx_neg,0,0,0],V);
        
        IO = [eye(n), zeros(n,mz)];
        OI = [zeros(mz,n), eye(mz)];
        
        Zy_neg=newlmi;
        lmiterm([Zy_neg,1,1,Y],[Ab;Cz],IO,'s');
        lmiterm([Zy_neg,1,1,Gamma_z],-OI',OI);
        lmiterm([Zy_neg,1,1,0],-OI'*Gamma0_z*OI);
        lmiterm([Zy_neg,1,2,0],[Bw; DDzw_fix]);
        if sum(strucDDzw_var(:))~=0; lmiterm([Zy_neg,1,2,DDzw_var],[zeros(n,mw); eye(mw)],1); end;
        lmiterm([Zy_neg,2,2,Gamma_w],-1,1);
        lmiterm([Zy_neg,2,2,0],-Gamma0_w);
        lmiterm([Zy_neg,0,0,0],W);
        
      
        Zxy_pos=newlmi;
        I_XY_pos=[eye(n-m) zeros(n-m,m)];
        lmiterm([-Zxy_pos 1 1 Y],1,1);
        lmiterm([-Zxy_pos 1 2 0],I_XY_pos');
        lmiterm([-Zxy_pos 2 2 X11],1,1);
        lmiterm([Zxy_pos  1 1 0],1e-3);
        
        LMIs=getlmis;
        
        % Step 1.4: Specifying the objective
        c = zeros(decnbr(LMIs),1);
        diagGamma_z = diag(decinfo(LMIs,Gamma_z));
        diagGamma_z = diagGamma_z(Iz1);
        c(diagGamma_z(diagGamma_z~=0)) = alpha((diagGamma_z~=0));
        
        % Step 1.5: Solving the SDP
        
        [fopt, xopt] = mincx(LMIs, c);
        
        LMIsopt = evallmi(LMIs, xopt);
        [Zx_neg, trash] = showlmi(LMIsopt, Zx_neg);
        [Zy_neg, trash] = showlmi(LMIsopt, Zy_neg);
        [trash, Zxy_pos] = showlmi(LMIsopt, Zxy_pos);
        check_LMIs = {'XY pos' , min(eig(Zxy_pos)) ;
            'X neg', max(eig(Zx_neg));
            'Y neg', max(eig(Zy_neg))}
        % Step 1.6: Retrieve variables
        X11 = dec2mat(LMIs, xopt,X11);
        X21 = dec2mat(LMIs, xopt,X21);
        X21(find(isnan(X21))) = eps;
        X = dec2mat(LMIs, xopt,X);
        Y = dec2mat(LMIs, xopt,Y);
        Gamma_z = dec2mat(LMIs, xopt,Gamma_z)+Gamma0_z;
        Gamma_w = dec2mat(LMIs, xopt,Gamma_w)+Gamma0_w;
        gamma=(diag(Gamma_z));
        XY_pos=[Y I_XY_pos';I_XY_pos X11];
        DDzw = DDzw_fix;
        if sum(strucDDzw_var(:))~=0; DDzw=DDzw+dec2mat(LMIs,xopt,DDzw_var); end;
        
        xopt(diagGamma_z(diagGamma_z~=0))=[];
end

% Step 1bis: Improve numerical conditioning
% -------------------------------------------
if options.numstab %Force eigenvalues of P further away from the imaginary axis    
    gamma = 1.01*gamma;
switch LMIparser
    
    case 'yalmip' %Parse with YALMIP, solve with SeDuMi
        
        % Step 1b.1: Declaring the variables
        X11 = sdpvar(n-m,n-m,'symmetric');
        X21= sdpvar(m, n-m);
        X= [X11; X21];
        Y = sdpvar(n,n,'symmetric');
        Gamma_w = sdpvar(mw);
        Gamma_z = sdpvar(mz);
        DDzw = sdpvar(mz,mw,'full');
        t=sdpvar(1);
        for i = 1:noc
            if alpha(i) == 0
                gamma(i) = 1;
            end
            
            Gamma_z(Iz1(i):Iz2(i),Iz1(i):Iz2(i)) = gamma(i) * eye(Mz(i)); % eye(Mz(i))
            Gamma_w(Iw1(i):Iw2(i),Iw1(i):Iw2(i)) = gamma(i) * eye(Mw(i)); % eye(Mw(i))
            DDzw(Iz1(i):Iz2(i),Iw1(i):Iw2(i)) = Dzw(Iz1(i):Iz2(i),Iw1(i):Iw2(i));
        end
        
        % Step 1b.2: Specifying the LMIs
       
        Zx = [Ab(:,1:n-m)'*X+X'*Ab(:,1:n-m), X'*Bw    , Cz(:,1:n-m)'   ;
            Bw'*X   , -Gamma_w, DDzw'  ;
            Cz(:,1:n-m)      , DDzw     , -Gamma_z];
        
        Zy = [Ab*Y+Y*Ab', Y*Cz' , Bw      ;
            Cz*Y    , -Gamma_z, DDzw     ;
            Bw'     , DDzw'  , -Gamma_w];
        
        I_XY_pos=[eye(n-m) zeros(n-m,m)];
        XY_pos=[Y I_XY_pos';I_XY_pos X11];
        constr = [XY_pos >= t , V'*Zx*V <= -1e-8 , W'*Zy*W <= -1e-8, t>=eps];
        
        
        % Step 1b.3: Specifying the objective
        goal = -t;
        
        
        % Step 1b.4: Solving the SDP
        sol = optimize(constr, goal, sdpsettings('solver',options.gammasolver))
        check_LMIs = {'P pos' , min(eig(double(XY_pos)));
            'Zx neg', max(eig(double(V'*Zx*V)))            ;
            'Zy neg', max(eig(double(W'*Zy*W)))            ;
            't', double(t)}
        
        % Step 1b.5: Retrieve variables
        X11 = double(X11);
        X21 = [double(X21)];
        X21(find(isnan(X21))) = eps;
        Y = double(Y);
        XY_pos=double(XY_pos);
        Gamma_z = double(Gamma_z);
        Gamma_w = double(Gamma_w);
        DDzw = double(DDzw);
        t=double(t);
        
    case 'lmilab' %Parse with LMILab, solve with LMILab
        
         setlmis([])
        
        % Step 1.2b : Declaring the variables
        [X11,~,sX11] = lmivar(1,[n-m,1]);
        [X21,~,sX21] = lmivar(2,[m,n-m]);
        [X,~,sX]=lmivar(3,[sX11;sX21]);
        [Y,mvar] = lmivar(1,[n,1]);
        Gamma0_z = zeros(mz);
        Gamma0_w = zeros(mw);
        strucGamma_w = zeros(mw);
        strucGamma_z = zeros(mz);
        DDzw_fix = zeros(mz,mw);
        strucDDzw_var=zeros(mz,mw);
        for i = 1:noc
            if alpha(i) == 0
                Gamma0_z(Iz1(i):Iz2(i),Iz1(i):Iz2(i)) = eye(Mz(i));
                Gamma0_w(Iw1(i):Iw2(i),Iw1(i):Iw2(i)) = eye(Mw(i));
            else
                Gamma0_z(Iz1(i):Iz2(i),Iz1(i):Iz2(i)) = gamma(i);
                Gamma0_w(Iw1(i):Iw2(i),Iw1(i):Iw2(i)) = gamma(i);
            end
            temp = reshape(mvar+(1:Mz(i)*sum(Mz(i+1:end)))', Mz(i), sum(Mz(i+1:end)));
            strucGamma_z(Iz1(i):Iz2(i),Iz2(i)+1:end) = temp;
            mvar = mvar + Mz(i)*sum(Mz(i+1:end));
            
            temp = reshape(mvar+(1:Mw(i)*sum(Mw(i+1:end)))', Mw(i), sum(Mw(i+1:end)));
            strucGamma_w(Iw1(i):Iw2(i),Iw2(i)+1:end) = temp;
            mvar = mvar + Mw(i)*sum(Mw(i+1:end));
            
            DDzw_fix(Iz1(i):Iz2(i),Iw1(i):Iw2(i)) = Dzw_hat(Iz1(i):Iz2(i),Iw1(i):Iw2(i));
            strucDDzw_var(Iz1(i):Iz2(i),Iw1(i):Iw2(i))=1;
        end
        
        [strucGamma_z]=strucGamma_z+strucGamma_z';
        [Gamma_z,~,sGz] = lmivar(3, strucGamma_z);
        
        [strucGamma_w]=strucGamma_w+strucGamma_w';
        [Gamma_w,~,sGw] = lmivar(3, strucGamma_w);
        
        strucDDzw_comp=find(strucDDzw_var);
        strucDDzw_var=zeros(mz,mw);
        strucDDzw_diff=setdiff([1:numel(strucDDzw_var)],strucDDzw_comp);
        strucDDzw_var(strucDDzw_diff)=[mvar+(1:numel(strucDDzw_diff))];
        
        [DDzw_var,~,sDDzwv]=lmivar(3,strucDDzw_var);
        
        [t,mt]=lmivar(1,[1,1]);
        
        % Step 1.3b : Specifying the LMIs
        
        IO=[eye(n-m) zeros(n-m,mw)];
        OI=[zeros(mw,n-m) eye(mw)];
        
        Zx_neg=newlmi;
        lmiterm([Zx_neg,1,1,X],[Ab(:,1:n-m) Bw]',IO,'s');
        lmiterm([Zx_neg,1,1,Gamma_w],-OI',OI);
        lmiterm([Zx_neg,1,1,0],-OI'*Gamma0_w*OI);
        lmiterm([Zx_neg,2,1,0],[Cz(:,1:n-m) DDzw_fix]);
        lmiterm([Zx_neg,2,1,DDzw_var],1,[zeros(n-m,mw); eye(mw)]');
        lmiterm([Zx_neg,2,2,Gamma_z],-1,1);
        lmiterm([Zx_neg,2,2,0],-Gamma0_z);
        lmiterm([Zx_neg,0,0,0],V);
        
        IO = [eye(n), zeros(n,mz)];
        OI = [zeros(mz,n), eye(mz)];
        
        Zy_neg=newlmi;
        lmiterm([Zy_neg,1,1,Y],[Ab;Cz],IO,'s');
        lmiterm([Zy_neg,1,1,Gamma_z],-OI',OI);
        lmiterm([Zy_neg,1,1,0],-OI'*Gamma0_z*OI);
        lmiterm([Zy_neg,1,2,0],[Bw; DDzw_fix]);
        lmiterm([Zy_neg,1,2,DDzw_var],[zeros(n,mw); eye(mw)],1);
        lmiterm([Zy_neg,2,2,Gamma_w],-1,1);
        lmiterm([Zy_neg,2,2,0],-Gamma0_w);
        lmiterm([Zy_neg,0,0,0],W);
        
        
        Zxy_pos=newlmi;
        I_XY_pos=[eye(n-m) zeros(n-m,m)];
        lmiterm([-Zxy_pos 1 1 Y],1,1);
        lmiterm([-Zxy_pos 1 2 0],I_XY_pos');
        lmiterm([-Zxy_pos 2 2 X11],1,1);
        lmiterm([Zxy_pos 1 1 t],1,1);
        
        t_pos=newlmi;
        lmiterm([-t_pos 1 1 t],1,1);
        lmiterm([-t_pos 2 2 t],1,1);
        
        LMIs=getlmis;
        
        % Step 1.4b : Specifying the objective
        c = zeros(decnbr(LMIs),1);
        c(mt) = -1;
        
        % Step 1.5b : Solving the SDP
        [fopt, xopt] = mincx(LMIs, c,[],xopt);
        
        LMIsopt = evallmi(LMIs, xopt);
        [Zx_neg, trash] = showlmi(LMIsopt, Zx_neg);
        [Zy_neg, trash] = showlmi(LMIsopt, Zy_neg);
        [trash, Zxy_pos] = showlmi(LMIsopt, Zxy_pos);
%        [trash, t_pos] = showlmi(LMIsopt, t_pos);

        check_LMIs = {'XY pos' , min(eig(Zxy_pos)) ;
            'X neg', max(eig(Zx_neg));
            'Y neg', max(eig(Zy_neg));
 %           't_pos', min(eig(t_pos))
 }
        
        % Step 1.6b : Retrieve variables
        X11 = dec2mat(LMIs, xopt,X11);
        X21 = dec2mat(LMIs, xopt,X21);
        X21(find(isnan(X21))) = eps;
        Y = dec2mat(LMIs, xopt,Y);
        Gamma_z = dec2mat(LMIs, xopt,Gamma_z)+Gamma0_z;
        Gamma_w = dec2mat(LMIs, xopt,Gamma_w)+Gamma0_w;
        XY_pos=[Y I_XY_pos';I_XY_pos X11];
        DDzw=dec2mat(LMIs,xopt,DDzw_var)+DDzw_fix;
        

end
end
% STEP 2: Reconstructing the controller
% -------------------------------------

%Step 2.1: Calculate M,N
[U,S,V] = svd(eye(n-m)-[X11 X21']*[Y(1:n-m,1:n-m);Y(n-m+1:n,1:n-m)]);
N1 = V*sqrt(S);
M = U*sqrt(S);
N2=-[Y(n-m+1:n,1:n-m) Y(n-m+1:end,n-m+1:end)]*[X11;X21]/(M');

%Step 3.2: Calculate Lyapunov matrix P trough reverse congruence transform
PI_c=[Y [eye(n-m);zeros(m,n-m)]; N1' N2' zeros(n-m)];
P=PI_c'\XY_pos/PI_c;

% Step 2.2: Retrieving the controller
% switch 1
%     case 1 %Controller parameters
nc = n-m;
Aa = blkdiag(A+beta*eye(n),beta*eye(nc));
Bwa = [Bw; zeros(nc,mw)];
Bua = [zeros(n,nc), Bu; eye(nc), zeros(nc,mu)];
Cza = [Cz, zeros(mz,nc)];
Dzua = [zeros(mz,nc), Dzu];
Cya = [zeros(nc,n), eye(nc); Cy, zeros(my,nc)];
Dywa = [zeros(nc,mw); Dyw];
Z0 = [Aa'*P+P*Aa, P*Bwa     , Cza';
    Bwa'*P      , -Gamma_w, DDzw'  ;
    Cza     , DDzw     , -Gamma_z];
L = [P*Bua; zeros(mw,nc+mu); Dzua]';
R = [Cya, Dywa, zeros(nc+my,mz)];

%     case 2  % -> transformed controller parameters. BUGGED, use exact controller parameters (case 1) instead
%         nc = n-m;
%         
%         X1=[X11 X21']
%         PI_I=[eye(n-m);zeros(m+no,n-m)]
%         
%         Z0 = [A_hat*Y+Y*A_hat', A_hat*PI_I       , Bw_hat      , Y*Cz_hat' ;
%             PI_I'*A_hat'      , X1*A_hat*PI_I+PI_I'*A_hat'*X1', X1*Bw_hat    , PI_I'*Cz_hat'   ;
%             Bw_hat'     , Bw_hat'*X1'   , -Gamma_w, DDzw'  ;
%             Cz_hat*Y    , Cz_hat*PI_I      , DDzw     , -Gamma_z];
%         L = [zeros(nc,n+no), eye(nc)     , zeros(nc,mw) , zeros(nc,mz);
%             Bu_hat'     , zeros(mu,nc), zeros(mu,mw), Dzu_hat'       ];
%         R = [eye(n+no)     , zeros(n+no,nc), zeros(n+no,mw), zeros(n+no,mz) ;
%             zeros(my,n+no), Cy_hat*PI_I      , Dyw_hat        , zeros(my,mz)];
% end

switch options.controllersolver
    case 'basiclmi' %Basiclmi
        theta = basiclmi(Z0,L,R,'Xmin');
    case 'lmilab' %LMIlab
        setlmis([]);
        theta=lmivar(2,[nc+mu nc+my]);
        theta_calc=newlmi;
        lmiterm([theta_calc 1 1 0],Z0);
        lmiterm([theta_calc 1 1 theta],L',R,'s')
        LMIs=getlmis;
        [fopt, xopt] = feasp(LMIs);
    
        theta = dec2mat(LMIs, xopt, theta);
        
    otherwise       % -> yalmip + SDPsolver
        theta = sdpvar(nc+mu,nc+my,'full');
        t = sdpvar(1);
        constr = [Z0+L'*theta*R+R'*theta'*L <= t*eye(n+nc+mw+mz)];
        options = sdpsettings('solver',options.controllersolver);
        goal = t;
        sol = solvesdp(constr, goal, options);
        theta = double(theta);
        check_t = value(t)
end
check_LMI = max(eig(Z0+L'*theta*R+R'*theta'*L))
Ac0 = theta(1:nc,1:nc);
Bc0= theta(1:nc,nc+1:end);
Cc0= theta(nc+1:end,1:nc);
Dc0 = theta(nc+1:end,nc+1:end);

if norm(Dyu,1) > 0,
        if norm(Dk0,1) > 0,
            Myuk = eye(dims.my)+Dyu*Dc0; 
            Mkyu = eye(dims.mu)+Dc0*Dyu;
            if svds(Myuk,1,'smallest') < 1e-6
                error('The controller could not be reconstructed due to a singularity. This is usually caused by an algebraic loop in the generalized plant.');
            else
                Ck = Mkyu\Cc0;
                Dk = Mkyu\Dc0;
                Ak = Ak0-Bk0*Dyu*Ck;
                Bk = Bk0/Myuk;
            end
        else
            Ak = Ac0-Bc0*Dyu*Cc0;
        end
    else
         Ak = Ac0; Bk = Bc0; Ck = Cc0; Dk = Dc0;
end

K=ss(Ac,Bc,Cc,Dc);
if ~isstable(lft(sys_P,K));
    warning('closed loop has not been stabilized. Controller design failed.')
end
end

function [ T,Cy_T,Dyw_T,m,index_r ] = meas_state_trans( Cy,Dyw )
%MEAS_STATE_TRANS Returns a state space transformation which transforms
%the measuremtent matrix Cy to the following structure Cy_T=[C1 C2;0 I]
%   
%   Inputs
%       Cy: measurement matrix
%       Dyw: Disturbance to meas matrix.
%
%   Output:
%       T: Corresponding state space transformation
%       Cy_T: Transformed measurement matrix
%       Dyw_T: Transformed disturbance to measurements matrix
%       m: Number of order reductions possible


%Pre-processing
[my,n]=size(Cy);
mw=length(Dyw(1,:));

% if rank(Cy)<my
%     error('Measurements are dependent')
% end

%Initialisation
T=eye(n);
row_T=[1:n]; %Stores available rows of T

%Transform Cy such that each row containts:
%Only zeros
%One '1'
%Leave the row the same when a disturbance is included in the measurement
%Results in matrix T1
for i=1:my
    if sum(abs(Cy(i,:)))~=0 && sum(abs(Dyw(i,:)))==0
        index=find(Cy(i,:));
        j=1;
        while ~ismember(index(j),row_T)
            j=j+1;
        end
        T(index(j),:)=Cy(i,:);
        row_T(index(j))=0;
    end
end
T1=T;

%Transform Cy such that the disturbance rows are on top
index=[];
for i=1:length(Dyw(:,1))
    if sum(abs(Dyw(i,:)))~=0
        index=[index;i];
    end
end
j=1;
for i=1:length(index)
    Cy([index(i) j],:)=Cy([j index(i)],:);
    Dyw([index(i) j],:)=Dyw([j index(i)],:);
    j=j+1;
end
Dyw_T=Dyw;
mi=length(index);
index_r=index;

%Transfrom Cy to Cy=[C1 C2;0 I]
j=n-(my-mi)+1;
Cy_T=Cy/(T);
T2=eye(n);
for i=mi+1:my
    
    index=find(Cy_T(i,:)<1+eps & Cy_T(i,:)>1-eps);
    
    T2([j index],:)=T2([index j],:);
    
    %T(:,[index j])=T(:,[j index]);
    j=j+1;
end
T=T2*T1;
Cy_T=Cy/(T);
m=my-mi;
end



