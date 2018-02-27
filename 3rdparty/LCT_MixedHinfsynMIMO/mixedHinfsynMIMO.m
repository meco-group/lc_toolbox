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

function [K, gamma, info] = mixedHinfsynMIMO(P,my,mu,alpha,beta,ch,opts)

    % mixedHinfsynMIMO - Mixed sensitivity H-infinity loop shaping tool for
    % solving multi-objective H-infinity norm minimization problems subject
    % to multiple constraints. 
    % 
    % // NOT YET BUT PLANNED: The function supports both unstable and
    % improper weighting functions and automatically applies free controller 
    % order reductions (i.e. without increasing the objective value). 
    % 
    % Laurens Jacobs, MECO Research Team, KU Leuven, July 2017
    % laurens.jacobs@kuleuven.be
    % 
    % This file was designed for integration in LCToolbox.
    % www.github.com/meco-group/LCToolbox 
    %
    % This function is inspired by earlier implementations of Michiel 
    % Dhadamus, Gijs Hilhorst and Goele Pipeleers.
    
    disp('##################################');
    disp('This is mixedHinfsynMIMO 1.1 beta.');
    disp('##################################');
    disp(' ');

    %% Preprocessing

    disp('Preprocessing...');
    
    % parse options
    if nargin < 6
        error('Not enough input arguments.');
    elseif nargin < 7
        opts = struct();
    end
    
    if ~isfield(opts,'gammasolver');
        opts.gammasolver = 'mosek';  
    end
    if ~isfield(opts,'controllersolver');
        opts.controllersolver = 'basiclmi';
    end
    if ~isfield(opts,'verbose');
        opts.verbose = 'off';
    end
    
    %  store dimensions
    [A,B,C,D] = ssdata(balreal(P)); % retrieve state-space matrices
    dims.n = length(A);             % plant's order
    dims.my = my;                   % number of measured outputs
    dims.mu = mu;                   % number of control inputs
    dims.noc = length(alpha == 0); 	% number of constraints
    assert(length(alpha) == length(ch), 'The number of weighting factors does not match the number of channels.');
    assert(all(alpha >= 0) && all(beta >= 0), 'Weighting factors can only have nonnegative values.');
    assert(isequal((alpha~=0),(beta==0)),'alpha and beta are inconsistent.'); 
    
    %  restructure the generalized plant 
    Bn = B(:,end-mu+1:end); Cn = C(end-my+1:end,:); Dn = D(end-my+1:end,end-mu+1:end);
    idxu = (size(D,2)-mu+1):1:size(D,2);
    idxy = (size(D,1)-my+1):1:size(D,1);
    for i = length(ch):-1:1
       dims.Mw(i) = length(ch(i).In);
       dims.Mz(i) = length(ch(i).Out);
       for j = flip(ch(i).In)
           Bn = [B(:,j) Bn];
           Dn = [D(idxy,j) Dn];
           idxu = [j idxu];
       end
       for j = flip(ch(i).Out)
           Cn = [C(j,:) ; Cn];
           Dn = [D(j,idxu) ; Dn];
           idxy = [j idxy];
       end
    end
    B = Bn; C = Cn; D = Dn;
    
    % reordered single inputs and outputs
    dims.mz = length(C(:,1))-dims.my;         % number of single output channels
    dims.mw = length(B(1,:))-dims.mu;         % number of single input channels
    
    % partition plant
    part.A = A;                       part.Bw = B(:,1:dims.mw);                  part.Bu = B(:,dims.mw+1:end);
    part.Cz = C(1:dims.mz,:);         part.Dzw = D(1:dims.mz,1:dims.mw);         part.Dzu = D(1:dims.mz,dims.mw+1:end);
    part.Cy = C(dims.mz+1:end,:);     part.Dyw = D(dims.mz+1:end,1:dims.mw);     part.Dyu = D(dims.mz+1:end,dims.mw+1:end);
    
    %% Synthesis problem
    %  Find the Lyapunov matrices and the optimal gamma values
    %  using the elimination/projection lemma and the Lyapunov shaping paradigm
   
   if strcmp(opts.gammasolver,'lmilab')
        synresults = synthesis_LMILAB(part,dims,alpha,beta);
   else % solve through YALMIP, e.g. with MOSEK, SeDuMi, SDPT3,...
        synresults = synthesis_YALMIP(part,dims,alpha,beta,opts.gammasolver);
   end

    %% Analysis problem (controller reconstruction)
    % Find the controller matrices (analytical solution)
    
    K = reconstruction(part,dims,synresults,opts);
    
    %% Postprocessing 

    info = '';
    gamma = synresults.gamma;
    
    if ~isstable(lft(P,K))
        warning('The closed loop could not be stabilized due to numerical troubles. The result may not be reliable. Try to rescale your plant appropriately and set appropriate solvers.');
    end
    
end

function synresults = synthesis_LMILAB(part,dims,alpha,beta)

    disp('Formulating synthesis problem...');
    
    % define variables
    setlmis([]);
    
        % reserve the first scalar optimization variables for gamma
        if sum(alpha) ~= 0 % to avoid error about outputs which are not assigned
            [gamma,~,struc_gamma] = lmivar(2,[sum(alpha ~= 0), 1]);
        end
    
        % construct Lyapunov submatrices
        X = lmivar(1,[dims.n 1]);          % upper left block inverse Lyapunov matrix
        [Y,xn] = lmivar(1,[dims.n 1]);     % upper left block Lyapunov matrix
        
    	% construct Gamma_W0
            struc_Gamma_W0 = zeros(dims.mw);
            idx_Gamma_W0_constr = zeros(dims.mw);
            
            % define its structure
            j = 1;
            var = 1;
            for i = 1:length(dims.Mw)
                if alpha(i) ~= 0;   % an objective: add gamma on the diagonal
                    struc_Gamma_W0(j:j+dims.Mw(i)-1,j:j+dims.Mw(i)-1) = var*eye(dims.Mw(i));
                    var = var + 1;
                else                % a constraint: remember this was a constraint for later
                    idx_Gamma_W0_constr(j:j+dims.Mw(i)-1,j:j+dims.Mw(i)-1) = beta(i)*eye(dims.Mw(i));
                end
                if i > 1            % add slack optimization variables to reduce conservatism
                    slack = reshape(1:dims.Mw(i)*sum(dims.Mw(1:i-1)),[sum(dims.Mw(1:i-1)), dims.Mw(i)])' + xn;
                    struc_Gamma_W0(j:j+dims.Mw(i)-1,1:sum(dims.Mw(1:i-1))) = slack;
                    struc_Gamma_W0(1:sum(dims.Mw(1:i-1)),j:j+dims.Mw(i)-1) = slack';
                    xn = xn + dims.Mw(i)*sum(dims.Mw(1:i-1));
                end
                j = j + dims.Mw(i);
            end
            
            % define as an LMI variable
            Gamma_W0 = lmivar(3,struc_Gamma_W0);
            
    	% construct Gamma_Z0
            struc_Gamma_Z0 = zeros(dims.mz);
            idx_Gamma_Z0_constr = zeros(dims.mz);
            
            % define its structure
            j = 1;
            var = 1;
            for i = 1:length(dims.Mz)
                if alpha(i) ~= 0;   % an objective: add gamma on the diagonal
                    struc_Gamma_Z0(j:j+dims.Mz(i)-1,j:j+dims.Mz(i)-1) = var*eye(dims.Mz(i));
                    var = var + 1;
                else                % a constraint: remember this was a constraint for later
                    idx_Gamma_Z0_constr(j:j+dims.Mz(i)-1,j:j+dims.Mz(i)-1) = beta(i)*eye(dims.Mz(i));
                end
                if i > 1            % add slack optimization variables to reduce conservatism
                    slack = reshape(1:dims.Mz(i)*sum(dims.Mz(1:i-1)),[sum(dims.Mz(1:i-1)), dims.Mz(i)])' + xn;
                    struc_Gamma_Z0(j:j+dims.Mz(i)-1,1:sum(dims.Mz(1:i-1))) = slack;
                    struc_Gamma_Z0(1:sum(dims.Mz(1:i-1)),j:j+dims.Mz(i)-1) = slack';
                    xn = xn + dims.Mz(i)*sum(dims.Mz(1:i-1));
                end
                j = j + dims.Mz(i);
            end
            
            % define as an LMI variable
            Gamma_Z0 = lmivar(3,struc_Gamma_Z0);
            
        % construct DL
            struc_DL = zeros(dims.mz,dims.mw);
            
            % define its structure
            j = dims.Mz(1)+1;
            for i = 2:length(dims.Mz)
                slack = reshape(1:dims.Mz(i)*sum(dims.Mw(1:i-1)),[sum(dims.Mw(1:i-1)), dims.Mz(i)])' + xn;
                struc_DL(j:j+dims.Mz(i)-1,1:sum(dims.Mw(1:i-1))) = slack;
                xn = xn + dims.Mz(i)*sum(dims.Mw(1:i-1));
                j = j + dims.Mw(i);
            end
       
            j = dims.Mw(1)+1;
            for i = 2:length(dims.Mw)
                slack = reshape(1:dims.Mw(i)*sum(dims.Mz(1:i-1)),[sum(dims.Mz(1:i-1)), dims.Mw(i)]) + xn;
                struc_DL(1:sum(dims.Mz(1:i-1)),j:j+dims.Mw(i)-1) = slack;
                xn = xn + dims.Mw(i)*sum(dims.Mz(1:i-1));
                j = j + dims.Mz(i);
            end
            
            % define as an LMI variable
            DLT = lmivar(3,struc_DL');
            DLconst = zeros(size(struc_DL));
            
            % partitioning of Dzw into the constant part of DL
            idcs = [];
            for i = 1:length(dims.Mw)
                active = ones(dims.Mz(i),dims.Mw(i));
                idcs = blkdiag(idcs, active);
            end
            idcs = find(idcs==1);
            DLconst(idcs) = part.Dzw(idcs);
            
    % calculate the null spaces for the projection/elimination lemma
        % only calculate the parts we know that won't be zero

        W12 = null([part.Cy  part.Dyw]);
        V12 = null([part.Bu' part.Dzu']);

        % cast into outer factors with appropriate dimensions

          W = blkdiag(W12, eye(dims.mz));
          V = [V12(1:dims.n, :) zeros(dims.n, dims.mw) ; zeros(dims.mw, size(V12,2)) eye(dims.mw) ; V12(dims.n+1:end, :) zeros(size(V12,1)-dims.n, dims.mw)];
          
    % elaborate the actual LMIs
        % first inequality: ZX < 0
       	
        ZX = newlmi;
        lmiterm([ZX 0 0 0],V);
        lmiterm([ZX 1 1 X],part.A,1,'s');
        lmiterm([ZX 1 2 0],part.Bw);
        lmiterm([ZX 1 3 X],1,part.Cz');
        if ~isempty(Gamma_W0); lmiterm([ZX 2 2 Gamma_W0],-1,1); end;
        lmiterm([ZX 2 2 0],-idx_Gamma_W0_constr);
        if ~isempty(DLT); lmiterm([ZX 2 3 DLT],1,1); end;
        lmiterm([ZX 2 3 0],DLconst');
        if ~isempty(Gamma_Z0); lmiterm([ZX 3 3 Gamma_Z0],-1,1); end;
        lmiterm([ZX 3 3 0],-idx_Gamma_Z0_constr);
               
        % second inequality: ZY < 0
        
        ZY = newlmi;
        lmiterm([ZY 0 0 0],W);
        lmiterm([ZY 1 1 Y],part.A',1,'s');
        lmiterm([ZY 1 2 Y],1,part.Bw);
        lmiterm([ZY 1 3 0],part.Cz');
        if ~isempty(Gamma_W0); lmiterm([ZY 2 2 Gamma_W0],-1,1); end;
        lmiterm([ZY 2 2 0],-idx_Gamma_W0_constr);
        if ~isempty(DLT); lmiterm([ZY 2 3 DLT],1,1); end;
        lmiterm([ZY 2 3 0],DLconst');
        if ~isempty(Gamma_Z0); lmiterm([ZY 3 3 Gamma_Z0],-1,1); end;
        lmiterm([ZY 3 3 0],-idx_Gamma_Z0_constr);
        
        % third inequality: P > 0
        
        PXY = newlmi;
        lmiterm([-PXY 1 1 X],1,1);
        lmiterm([-PXY 2 2 Y],1,1);
        lmiterm([-PXY 1 2 0],1);
        
    % formulate the objective
    LMIs = getlmis;                             % get the LMI system
    c = zeros(decnbr(LMIs),1);                  % get the variables of the problem
    c(1:sum(alpha~=0)) = alpha(alpha ~= 0);     % objective = c'*x

    % solve
    disp('Solving synthesis problem...');
    T = evalc('[~, xopt] = mincx(LMIs, c);');   % evalc to hide ugly output
    if strfind(T,'INFEASIBLE'); error('The optimization failed: the constraints were found infeasible.'); end;
    
    % extract variables 
    synresults.X = dec2mat(LMIs,xopt,X);
    synresults.Y = dec2mat(LMIs,xopt,Y);
    
    synresults.gamma = [];
    synresults.Gamma_W0 = idx_Gamma_W0_constr;
    synresults.Gamma_Z0 = idx_Gamma_Z0_constr;
    synresults.DL = DLconst;
    
    if sum(struc_gamma(:))~=0
        synresults.gamma = dec2mat(LMIs,xopt,gamma); 
    end
    if sum(struc_Gamma_W0(:))~=0
        synresults.Gamma_W0 = dec2mat(LMIs,xopt,Gamma_W0) + idx_Gamma_W0_constr;
    end
    if sum(struc_Gamma_Z0(:))~=0
        synresults.Gamma_Z0 = dec2mat(LMIs,xopt,Gamma_Z0) + idx_Gamma_Z0_constr;  
    end
    if sum(struc_DL(:))~=0
        synresults.DL = dec2mat(LMIs,xopt,DLT)' + DLconst;
    end
    
end

function synresults = synthesis_YALMIP(part,dims,alpha,beta,solver)

    disp('Formulating synthesis problem...');

    % define variables 
    
    X = sdpvar(dims.n);
    Y = sdpvar(dims.n);
    gamma = sdpvar(1,sum(alpha~=0));
    
    % construct matrices for LMI formulation
    
    
    idx_ga = cumsum(alpha~=0);
    
    Gamma_W0 = sdpvar(dims.mw);
    for i = 1:length(dims.Mw)
        if alpha(i) ~= 0
            Gamma_W0(sum(dims.Mw(1:i-1))+(1:dims.Mw(i)),sum(dims.Mw(1:i-1))+(1:dims.Mw(i))) = gamma(idx_ga(i))*eye(dims.Mw(i));
        else 
            Gamma_W0(sum(dims.Mw(1:i-1))+(1:dims.Mw(i)),sum(dims.Mw(1:i-1))+(1:dims.Mw(i))) = beta(i)*eye(dims.Mw(i));
        end
    end
    Gamma_Z0 = sdpvar(dims.mz);
    for i = 1:length(dims.Mz)
        if alpha(i) ~= 0
            Gamma_Z0(sum(dims.Mz(1:i-1))+(1:dims.Mz(i)),sum(dims.Mz(1:i-1))+(1:dims.Mz(i))) = gamma(idx_ga(i))*eye(dims.Mz(i));
        else
            Gamma_Z0(sum(dims.Mz(1:i-1))+(1:dims.Mz(i)),sum(dims.Mz(1:i-1))+(1:dims.Mz(i))) = beta(i)*eye(dims.Mz(i));
        end
    end
    DL = sdpvar(sum(dims.Mz),sum(dims.Mw),'full');
    for i = 1:length(dims.Mw)
            DL(sum(dims.Mz(1:i-1))+(1:dims.Mz(i)),sum(dims.Mw(1:i-1))+(1:dims.Mw(i))) = part.Dzw(sum(dims.Mz(1:i-1))+(1:dims.Mz(i)),sum(dims.Mw(1:i-1))+(1:dims.Mw(i)));
    end
    
    % calculate null spaces
    
      W12 = null([part.Cy  part.Dyw]);
      V12 = null([part.Bu' part.Dzu']);
      W = blkdiag(W12, eye(dims.mz));
      V = [V12(1:dims.n, :) zeros(dims.n, dims.mw) ; zeros(dims.mw, size(V12,2)) eye(dims.mw) ; V12(dims.n+1:end, :) zeros(size(V12,1)-dims.n, dims.mw)];
          
    % elaborate the actual LMIs
    
    Zx = [  part.A*X+X*part.A'        part.Bw     X*part.Cz'   ;
            part.Bw'                  -Gamma_W0   DL'          ;
            part.Cz*X                 DL          -Gamma_Z0    ];
      
    Zy = [  part.A'*Y+Y*part.A        Y*part.Bw   part.Cz'     ;
            part.Bw'*Y                -Gamma_W0   DL'          ;
            part.Cz                   DL          -Gamma_Z0    ];
        
    ZX = V'*Zx*V;
    ZY = W'*Zy*W;
    PXY = [X eye(dims.n) ; eye(dims.n) Y];
    constr = (PXY >= 0) + (ZX <= 0) + (ZY <= 0);
    
    % formulate the objective
    
    c = alpha(alpha~=0);
    obj = gamma*c(:);
    
    % solve the problem
    
    disp('Solving synthesis problem...');
    
    options = sdpsettings('solver',solver);
    sol = optimize(constr, obj, options);
    if strfind(sol.info,'Infeasible'); error('The optimization failed: the constraints were found infeasible.'); end;
    
    % extract the variables
    
    synresults.X = double(X);
    synresults.Y = double(Y);
    synresults.gamma = double(gamma);
    synresults.Gamma_W0 = double(Gamma_W0);
    synresults.Gamma_Z0 = double(Gamma_Z0);
    synresults.DL = double(DL);
    
end

function K = reconstruction(part,dims,synresults,opts)

    disp('Reconstructing the controller...');
    
    % reconstruct the the Lyapunov matrix  
    
    [U,S] = svd(eye(dims.n)-synresults.X*synresults.Y);
    M = U*sqrt(S);
    PI = [synresults.X eye(dims.n) ; M' zeros(dims.n)];
    Ptransf = [synresults.X eye(dims.n) ; eye(dims.n) synresults.Y];
    P = PI'\Ptransf/PI;
    
    % relax gamma values with 1% to ease the reconstruction

    relax = 0.01;
    Gamma_W0 = synresults.Gamma_W0 + relax*diag(diag(synresults.Gamma_W0));
    Gamma_Z0 = synresults.Gamma_Z0 + relax*diag(diag(synresults.Gamma_Z0));
    
    % substitute the Lyapunov matrix parts in the BMI formulation and cast 
    % into this format: Z + U'*theta*V + V'*theta'*U < 0
    
    A0 = zeros(2*dims.n,2*dims.n); A0(1:dims.n,1:dims.n) = part.A;
    B0 = [part.Bw ; zeros(dims.n,dims.mw)];
    C0 = [part.Cz zeros(dims.mz,dims.n)];
    D0 = synresults.DL;

    Z = [A0'*P+P*A0   P*B0          C0'       ;
         B0'*P        -Gamma_W0     D0'       ; 
         C0           D0            -Gamma_Z0 ];  

    Buhat = [zeros(dims.n) part.Bu ; eye(dims.n) zeros(dims.n, dims.mu)];
    Dzuhat = [zeros(dims.mz,dims.n) part.Dzu];
    Cyhat = [zeros(dims.n) eye(dims.n) ; part.Cy zeros(dims.my,dims.n)];
    Dyhat = [zeros(dims.n,dims.mw) ; part.Dyw];

    U = [P*Buhat ; zeros(length(synresults.Gamma_W0),size(Dzuhat,2)) ; Dzuhat]';
    V = [Cyhat Dyhat zeros(size(Cyhat,1),size(synresults.Gamma_Z0,1))];
    
    % solve for the controller parameters
    
    switch opts.controllersolver
        case 'basiclmi'
            theta = basiclmi(Z,U,V,'Xmin,Shift');
        case 'lmilab'
            setlmis([]);
            theta = lmivar(2,[dims.n+dims.mu dims.n+dims.my]);
            receq = newlmi;
            lmiterm([receq 1 1 0],Z);
            lmiterm([receq 1 1 theta],U',V,'s')
            LMIs = getlmis;
            T = evalc('[~, xopt] = feasp(LMIs);');  % evalc to hide ugly output
            if strfind(T,'INFEASIBLE'); error('The reconstruction problem appears to be infeasible. The synthesis might have returned ill-conditioned results. Try to use ''basiclmi'', which is a numerically more stable reconstruction algorithm.'); end;
            theta = dec2mat(LMIs,xopt,theta);
        otherwise % through YALMIP
            theta = sdpvar(dims.n+dims.mu, dims.n+dims.my,'full');
            constr = [Z + U'*theta*V + V'*theta'*U < 0];
            options = sdpsettings('solver',opts.controllersolver);
            optimize(constr,[],options);
            theta = double(theta);
    end
    
    Ak0 = theta(1:dims.n,1:dims.n);
    Bk0 = theta(1:dims.n,dims.n+1:end);
    Ck0 = theta(dims.n+1:end,1:dims.n);
    Dk0 = theta(dims.n+1:end,dims.n+1:end);
    
    % compensate for feedthrough matrix Dyu

    if norm(part.Dyu,1) > 0,
        if norm(Dk0,1) > 0,
            Myuk = eye(dims.my)+part.Dyu*Dk0; 
            Mkyu = eye(dims.mu)+Dk0*part.Dyu;
            if svds(Myuk,1,'smallest') < 1e-6
                error('The controller could not be reconstructed due to a singularity. This is usually caused by an algebraic loop in the generalized plant.');
            else
                Ck = Mkyu\Ck0;
                Dk = Mkyu\Dk0;
                Ak = Ak0-Bk0*part.Dyu*Ck;
                Bk = Bk0/Myuk;
            end
        else
            Ak = Ak0-Bk0*part.Dyu*Ck0;
            Bk = Bk0; Ck = Ck0; Dk = Dk0;
        end
    else
         Ak = Ak0; Bk = Bk0; Ck = Ck0; Dk = Dk0;
    end
    
    % return the controller
    
    K = ss(Ak,Bk,Ck,Dk);

end