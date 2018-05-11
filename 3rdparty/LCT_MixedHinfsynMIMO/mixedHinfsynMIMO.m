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

function [K, gamma, info] = mixedHinfsynMIMO(P,Wi,Wo,my,mu,alpha,beta,ch,opts)

    % mixedHinfsynMIMO - Mixed sensitivity controller synthesis tool for
    % solving multi-objective H-infinity norm minimization problems subject
    % to multiple H-infinity norm constraints.
    %
    % Laurens Jacobs, MECO Research Team, KU Leuven, May 2018
    % laurens.jacobs@kuleuven.be
    % 
    % Inputs: 
    %   P:      generalized plant 
    %   Wi:     unstable and/or impulsive input filter (only unstable and/or impulsive modes are allowed, include other modes in P)
    %   Wo:     unstable and/or impulsive output filter (only unstable and/or impulsive modes are allowed, include other modes in P)
    %   my:     number of measurements 
    %   mu:     number of control inputs
    %   alpha:  weighting factor for the norms in the objective (if a constraint: 0)
    %   beta:   constraint on the norm (if not a constraint: 0) 
    %   ch:     defines the inputs and outputs of every performance channel with indices,
    %           e.g. ch(1).In = [1 2], ch(1).Out = [1]
    %   opts:   allows to change default solver settings
    %
    % Outputs:
    %   K:      optimal controller, which is the result of the following
    %           optimization problem:
    %
    %                        +----+     +-----+     +----+
    %                 w  +-->+ Wi +---->+     +---->+ Wo +--->  z
    %                        +----+     |  P  |     +----+
    %                               +-->+     +---+
    %                               |   +-----+   |
    %                            u  |             |  y
    %                               |   +-----+   |
    %                               +---+  K  +<--+
    %                                   +-----+
    %
    %            minimize (alpha(1)*gamma(1) + alpha(2)*gamma(2) + ...)
    %               K
    %
    %            subject to |z(ch(1).Out)/w(ch(1).In)| <= (alpha(1)~=0)*gamma(1) + beta(1) 
    %                       |z(ch(2).Out)/w(ch(2).In)| <= (alpha(2)~=0)*gamma(2) + beta(2) 
    %                       ...
    %
    %   gamma:  gamma(:,1) are the anticipated optimal gamma (if objective term) or beta (if constraint) values (synthesis results)
    %           gamma(:,2) are the actual achieved gamma (if objective term) or beta (if constraint) values (reconstruction results)
    % 
    % References:
    %   - P. Gahinet and P. Apkarian. "A Linear Matrix Inequality Approach
    %     to H-infinity control." International Journal of Robust Nonlinear
    %     Control, vol. 4, 1994.
    %   - H. Koroglu. "H-infinity Synthesis with Unstable Weighting
    %     Filters: An LMI Solution." 52nd IEEE Conference on Decision and
    %     Control (CDC), Florence, Italy, 2013.
    %   - P. Gahinet, C. Scherer and M. Chilali. "Multiobjective
    %     Output-Feedback Control via LMI Optimization." IEEE Transactions 
    %     on Automatic Control, vol. 42, no. 7, 1997.
    %   - X. Xin, L. Guo and C. Feng. "Reduced-Order Controllers for
    %     Continuous and Discrete-time Singular H-infinity Control Problems
    %     Based on LMI." Automatica, vol. 32, no. 11, 1996.
    %   - T. Iwasaki and R.E. Skelton. "All Controllers for the General
    %     H-infinity Control Problem: LMI Existence Conditions and State
    %     Space Formulas." Automatica, vol. 30, no. 8, 1994.
    %   - Y. Feng and S. Ziyi. "H-infinity Control with Output Weights For
    %     Descriptor Systems: An LMI Approach." Proceedings of the 33rd
    %     Chinese Control Conference, Nanjing, China, 2014.
    %   - I. Masubuchi. "Output feedback controller synthesis for
    %     descriptor systems satisfying closed-loop dissipativity."
    %     Automatica, vol. 43, no. 2, 2007.

    %% Preprocessing
    % Restructure the generalized plant and check and store dimensions

    msgout(' ___________________________________\n| <strong>This is mixedHinfsynMIMO, v1.1.</strong>   |\n| Laurens Jacobs, May 2018        |\n|___________________________________|\n\n');
    msgout('<strong>Preprocessing...</strong>\n');
    
    % parse options
    if nargin < 8
        error('Not enough input arguments.');
    elseif nargin < 9
        opts = struct();
    end
    
    if ~isfield(opts,'gammasolver') || ~isfield(opts.gammasolver,'solver') || ~ischar(opts.gammasolver.solver)
        opts.gammasolver.solver = 'mosek';
    end
    
    if ~isfield(opts,'controllersolver') || ~isfield(opts.controllersolver,'solver')  || ~ischar(opts.controllersolver.solver)
        opts.controllersolver.solver = 'basiclmi';
    end
    
    % preprocess the generalized plant
    assert(size(Wo,2)+my == size(P,1) && size(Wi,1)+mu == size(P,2),'Unstable input and/or output weighting filter dimensions are not consistent.');
    assert(P.Ts == Wi.Ts && P.Ts == Wo.Ts && P.Ts >= 0, 'The sampling time of the unstable input and/or output weighting filter is different from the sampling time of the generalized plant.');
    [A,B,C,D] = ssdata(balreal(P));          % retrieve state-space matrices of stabilizible generalized plant
    [Ai,Bi,Ci,Di] = ssdata(balreal(Wi));     % retrieve state-space matrices of unstable input filter
    [Ao,Bo,Co,Do] = ssdata(balreal(Wo));     % retrieve state-space matrices of unstable output filter
    
    %  store dimensions
    stabplant.dims.n = length(A);            % stabilizable generalized plant's order
    usfilter.dims.ni = length(Ai);           % unstable input filter's order 
    usfilter.dims.no = length(Ao);           % unstable output filter's order
    stabplant.dims.my = my;                  % number of measured outputs
    stabplant.dims.mu = mu;                  % number of control inputs
    stabplant.isdiscrete = (P.Ts~=0);        % discrete time flag
    assert(length(alpha) == length(ch), 'The number of weighting factors does not match the number of channels.');
    assert(all(alpha >= 0) && all(beta >= 0), 'Weighting factors can only have nonnegative values.');
    assert(isequal((alpha~=0),(beta==0)),'alpha and beta are inconsistent.'); 
    if ~(isempty(pole(Wi)) || isempty(pole(Wo))); warning('You have both an unstable input filter and an unstable output filter. Redundant dynamics will therefore appear in the controller.'); end
    if stabplant.isdiscrete
        assert(all(abs(pole(Wi))>=1) && all(abs(pole(Wo))>=1), 'The unstable input and/or output weighting filter contains stable modes. All stable modes are supposed to be included in the generalized plant P.');
    else
        assert(all(real(pole(Wi))>=0) && all(real(pole(Wo))>=0), 'The unstable input and/or output weighting filter contains stable modes. All stable modes are supposed to be included in the generalized plant P.');
    end
    
    %  restructure the weighting filters
    Bni = zeros(size(Bi,1),0); Dni = zeros(size(Di,1),0);
    Cno = zeros(0,size(Co,2)); Dno = zeros(0,size(Do,2));
    for i = length(ch):-1:1
       usfilter.dims.Mw(i) = length(ch(i).In);
       usfilter.dims.Mz(i) = length(ch(i).Out);
       for j = flip(ch(i).In)
           % input filter
           Bni = [Bi(:,j) Bni];
           Dni = [Di(:,j) Dni];
       end
       for j = flip(ch(i).Out)
           % output filter
           Cno = [Co(j,:) ; Cno];
           Dno = [Do(j,:) ; Dno];
       end
    end
    
    % reordered single inputs and outputs
    usfilter.dims.mz = size(Cno,1);                       % number of single output channels
    usfilter.dims.mw = size(Bni,2);                       % number of single input channels
    stabplant.dims.mz = size(C,1)-stabplant.dims.my;      % number of signals between stabilizable plant and unstable output filter
    stabplant.dims.mw = size(B,2)-stabplant.dims.mu;      % number of signals between unstable input filter and stabilizable plant
    
    % partition the stabilizable generalized plant
    stabplant.A = A;                                 stabplant.Bw = B(:,1:stabplant.dims.mw);                            stabplant.Bu = B(:,stabplant.dims.mw+1:end);
    stabplant.Cz = C(1:stabplant.dims.mz,:);         stabplant.Dzw = D(1:stabplant.dims.mz,1:stabplant.dims.mw);         stabplant.Dzu = D(1:stabplant.dims.mz,stabplant.dims.mw+1:end);
    stabplant.Cy = C(stabplant.dims.mz+1:end,:);     stabplant.Dyw = D(stabplant.dims.mz+1:end,1:stabplant.dims.mw);     stabplant.Dyu = D(stabplant.dims.mz+1:end,stabplant.dims.mw+1:end);
    
    % gather state-space matrices of unstable weighting filters
    usfilter.Ai = Ai; usfilter.Bi = Bni; usfilter.Ci = Ci; usfilter.Di = Dni;
    usfilter.Ao = Ao; usfilter.Bo = Bo; usfilter.Co = Cno; usfilter.Do = Dno;
    
    %% Synthesis problem
    %  Find the Lyapunov matrices and the optimal gamma values (= the
    %  optimal closed loop) using the elimination/projection lemma and the 
    %  Lyapunov shaping paradigm
   
    [synresults, extplant, aux] = synthesis(stabplant,usfilter,alpha,beta,opts.gammasolver);

    %% Controller reconstruction
    % Reconstruct the values of the state-space matrices of the controller
    % for the optimal closed loop found in the synthesis problem
    
    if extplant.isdiscrete
        K = reconstruction_dt(stabplant,usfilter,extplant,synresults,aux,opts.controllersolver);
        K.Ts = P.Ts;
    else
        K = reconstruction_ct(stabplant,usfilter,extplant,synresults,aux,opts.controllersolver);
    end
    
    %% Analysis problem 
    % Check whether the reconstruction went fine by solving the closed-loop
    % analysis problem and compare the norms with the synthesis results
    
    spcl = lft(P,K);
    if ~isstable(spcl)
        warning('off','backtrace');
        warning('The reconstructed controller does not stabilize the generalized plant. Your problem is likely to be ill-conditioned and the results may thus be inaccurate.');
        warning('on','backtrace');
    end
    
    gamma = ones(length(alpha),2);
    gamma(alpha~=0,1) = synresults.gamma';
    gamma(beta~=0,1) = beta(beta~=0)';
    pcl = minreal(Wo*spcl*Wi,[],false);
    t = false; s = false;
    for i = 1:length(ch)
        gamma(i,2) = hinfnorm(pcl(ch(i).Out,ch(i).In),0.005);
        if beta(i)~=0 && gamma(i,2) > gamma(i,1)
            t = t | true; 
        end
        if alpha(i)~=1 && gamma(i,2)/gamma(i,1) > 1.2
            s = s | true;
        end
    end
    warning('off','backtrace');
    if t; warning('Although the synthesis problem seems feasible, the solution violates the constraints after the controller reconstruction. Your problem is likely to be ill-conditioned and the results may thus be inaccurate.'); end
    if s; warning('The performance achieved by the reconstructed controller differs more than 20% with respect to the predicted performance. Your problem is likely to be ill-conditioned and the results may thus be inaccurate.'); end
    warning('on','backtrace');
    
    %% Postprocessing 
    % No additional info is currently fed to the output
    
    info = ''; 
    msgout('<strong>mixedHinfsynMIMO was able to finish properly!</strong>\n\n');
    
end

function [synresults, extplant, aux] = synthesis(stabplant,usfilter,alpha,beta,gammasolver)

    msgout('<strong>Extending the generalized plant for including unstable weighting filters...</strong>\n');

    % construct auxiliary matrices for LMI formulation
    Ahat = [usfilter.Ao                                  usfilter.Bo*stabplant.Cz;
            zeros(stabplant.dims.n, usfilter.dims.no)    stabplant.A             ];
            
    Atilde = [stabplant.A                                  stabplant.Bw*usfilter.Ci;
              zeros(usfilter.dims.ni, stabplant.dims.n)    usfilter.Ai             ];
          
    Bwhat = [usfilter.Bo*stabplant.Dzw;
             stabplant.Bw             ];
    
    Bwtilde = [stabplant.Bw*usfilter.Di;
               usfilter.Bi             ];
    
    Czhat = [usfilter.Co usfilter.Do*stabplant.Cz];
    
    Cztilde = [stabplant.Cz stabplant.Dzw*usfilter.Ci];
    
    Buhat = [usfilter.Bo*stabplant.Dzu;
             stabplant.Bu             ];
    
    Cytilde = [stabplant.Cy stabplant.Dyw*usfilter.Ci];
    
    At = [Ahat Buhat ; Czhat usfilter.Do*stabplant.Dzu];                                                     % Sylvester-like equation: 
    Bt = -usfilter.Ai;                                                                                       % At*[PIihat;Gammai] + [PIihat;0]*Bt = Ct
    Ct = [Bwhat*usfilter.Ci ; usfilter.Do*stabplant.Dzw*usfilter.Ci];
    left = kron(eye(usfilter.dims.ni),At) + kron(Bt,blkdiag(eye(stabplant.dims.n+usfilter.dims.no),zeros(usfilter.dims.mw,stabplant.dims.mu)));
    right = Ct(:);
    if ~isempty(left); vars = lsqlin(left,right,[],[]); else; vars = zeros(size(left,1),size(right,2)); end
    vars = reshape(vars,[stabplant.dims.n+usfilter.dims.no+stabplant.dims.mu,usfilter.dims.ni]);                                                   
    PIihat = vars(1:stabplant.dims.n+usfilter.dims.no,:); Gammai = vars(end-stabplant.dims.mu+1:end,:);
    residual = At*[PIihat ; Gammai] + [PIihat ; zeros(usfilter.dims.mz,usfilter.dims.ni)]*Bt - Ct;          % check solution (avoid continuing while infeasible or ill-conditioned)
    assert(isempty(residual) || max(abs(reshape(residual,[numel(residual),1]))) < sqrt(eps), 'Could not find a suitable partitioning to eliminate the unstable input filter dynamics.');
    PIihat(abs(PIihat)<sqrt(eps)) = 0;
    
    At = [Atilde Bwtilde ; Cytilde stabplant.Dyw*usfilter.Di]';                                              % Sylvester-like equation: 
    Bt = -usfilter.Ao';                                                                                      % At*[PIotilde;Gammao] + [PIotilde;0]*Bt = Ct
    Ct = [usfilter.Bo*Cztilde usfilter.Bo*stabplant.Dzw*usfilter.Di]';
    left = kron(eye(usfilter.dims.no),At) + kron(Bt,blkdiag(eye(stabplant.dims.n+usfilter.dims.ni),zeros(usfilter.dims.mz,stabplant.dims.my)));
    right = Ct(:);
    if ~isempty(left); vars = lsqlin(left,right,[],[]); else; vars = zeros(size(left,1),size(right,2)); end
    vars = reshape(vars,[stabplant.dims.n+usfilter.dims.ni+stabplant.dims.my,usfilter.dims.no]);  
    PIotilde = vars(1:stabplant.dims.n+usfilter.dims.ni,:)'; Gammao = vars(end-stabplant.dims.my+1:end,:)';
    residual = At*[PIotilde' ; Gammao'] + [PIotilde' ; zeros(usfilter.dims.no,usfilter.dims.mw)']*Bt - Ct;  % check solution (avoid continuing while infeasible or ill-conditioned)
    assert(isempty(residual) || max(abs(reshape(residual,[numel(residual),1]))) < sqrt(eps), 'Could not find a suitable partitioning to eliminate the unstable output filter dynamics.');
    PIotilde(abs(PIotilde)<sqrt(eps)) = 0;
    
    PIo = PIotilde(:,1:stabplant.dims.n);
    THETAo = PIotilde(:,stabplant.dims.n+1:end);
    PIi = PIihat(end-stabplant.dims.n+1:end,:);
    THETAi = PIihat(1:end-stabplant.dims.n,:);
    Zohat = [(THETAo+THETAi-PIo*PIi) (PIo*stabplant.Bu-usfilter.Bo*stabplant.Dzu)];
    Zi = stabplant.Cy*PIi - stabplant.Dyw*usfilter.Ci;
    
    Ei = [zeros(stabplant.dims.mu,usfilter.dims.ni) eye(stabplant.dims.mu)];
    Eo = [zeros(stabplant.dims.my,usfilter.dims.no) eye(stabplant.dims.my)];
    
    extplant.A = [usfilter.Ao                                                     usfilter.Bo*Cztilde;
                  zeros(stabplant.dims.n+usfilter.dims.ni,usfilter.dims.no)       Atilde             ];
    extplant.Bw = [usfilter.Bo*stabplant.Dzw*usfilter.Di;
                   Bwtilde                              ];
    extplant.Bu = [-PIihat                Buhat                                    ;
                   eye(usfilter.dims.ni)  zeros(usfilter.dims.ni,stabplant.dims.mu)];
    extplant.Cz = [usfilter.Co usfilter.Do*Cztilde];
    extplant.Dzw = usfilter.Do*stabplant.Dzw*usfilter.Di;
    extplant.Dzu = usfilter.Do*stabplant.Dzu*Ei;
    extplant.Cy = [eye(usfilter.dims.no)                       -PIotilde; 
                   zeros(stabplant.dims.my,usfilter.dims.no)   Cytilde  ];
    extplant.Dyw = Eo'*stabplant.Dyw*usfilter.Di;
    
    extplant.dims.n = stabplant.dims.n+usfilter.dims.ni+usfilter.dims.no;
    extplant.dims.mw = usfilter.dims.mw;
    extplant.dims.mz = usfilter.dims.mz;
    extplant.dims.mu = stabplant.dims.mu + usfilter.dims.ni;
    extplant.dims.my = stabplant.dims.my + usfilter.dims.no;
    
    extplant.isdiscrete = stabplant.isdiscrete; 
    
    msgout('<strong>Checking for controller order reduction possibilities...</strong>\n');
    
    % detect free controller order reductions
    [extplant, nxy, nxu, Ty, Tu] = detectFCOR(extplant);
    A11X = extplant.A(1:end-nxy,1:end-nxy);
    A11Y = extplant.A(1:end-nxu,1:end-nxu);
    A12 = extplant.A(1:end-nxu,end-nxu+1:end);
    A21 = extplant.A(end-nxy+1:end,1:extplant.dims.n-nxy);
    Cz1 = extplant.Cz(:,1:end-nxy);
    Bw1 = extplant.Bw(1:end-nxu,:);

    msgout('<strong>Parsing synthesis problem into an SDP...</strong>\n');

    % define SDP variables 
    X11 = sdpvar(stabplant.dims.n+usfilter.dims.ni+usfilter.dims.no-nxy);           % part of the Lyapunov matrix
    X12 = sdpvar(stabplant.dims.n+usfilter.dims.ni+usfilter.dims.no-nxy,nxy);       % part of the Lyapunov matrix
    X22 = sdpvar(nxy,nxy);                                                          % part of the Lyapunov matrix
    X = [X11 X12 ; X12' X22];                                                       % for discrete-time we need the full block
    Y11 = sdpvar(stabplant.dims.n+usfilter.dims.no+usfilter.dims.ni-nxu);           % part of the inverse Lyapunov matrix
    Y21 = sdpvar(nxu,stabplant.dims.n+usfilter.dims.ni+usfilter.dims.no-nxu);       % part of the inverse Lyapunov matrix
    Y22 = sdpvar(nxu,nxu);                                                          % part of the inverse Lyapunov matrix
    Y = [Y11 Y21' ; Y21 Y22];                                                       % for discrete-time we need the full block
    gamma = sdpvar(1,sum(alpha~=0));                                                % H-infinity norms involved in the optimization problem
    
    idx_obj = cumsum(alpha~=0);
    
    Gamma_W = sdpvar(usfilter.dims.mw);
    for i = 1:length(usfilter.dims.Mw)
        if isempty(gamma)
            Gamma_W(sum(usfilter.dims.Mw(1:i-1))+(1:usfilter.dims.Mw(i)),sum(usfilter.dims.Mw(1:i-1))+(1:usfilter.dims.Mw(i))) = beta(i)*eye(usfilter.dims.Mw(i));
        else
            Gamma_W(sum(usfilter.dims.Mw(1:i-1))+(1:usfilter.dims.Mw(i)),sum(usfilter.dims.Mw(1:i-1))+(1:usfilter.dims.Mw(i))) = (alpha(i)*gamma(idx_obj(i))+beta(i))*eye(usfilter.dims.Mw(i));
        end
    end
    
    Gamma_Z = sdpvar(usfilter.dims.mz);
    for i = 1:length(usfilter.dims.Mz)
        if isempty(gamma)
            Gamma_Z(sum(usfilter.dims.Mz(1:i-1))+(1:usfilter.dims.Mz(i)),sum(usfilter.dims.Mz(1:i-1))+(1:usfilter.dims.Mz(i))) = beta(i)*eye(usfilter.dims.Mz(i));
        else
            Gamma_Z(sum(usfilter.dims.Mz(1:i-1))+(1:usfilter.dims.Mz(i)),sum(usfilter.dims.Mz(1:i-1))+(1:usfilter.dims.Mz(i))) = (alpha(i)*gamma(idx_obj(i))+beta(i))*eye(usfilter.dims.Mz(i));
        end
    end
    
    Dstack = sdpvar(sum(usfilter.dims.Mz),sum(usfilter.dims.Mw),'full');
    for i = 1:length(usfilter.dims.Mw)
        Dstack(sum(usfilter.dims.Mz(1:i-1))+(1:usfilter.dims.Mz(i)),sum(usfilter.dims.Mw(1:i-1))+(1:usfilter.dims.Mw(i))) = extplant.Dzw(sum(usfilter.dims.Mz(1:i-1))+(1:usfilter.dims.Mz(i)),sum(usfilter.dims.Mw(1:i-1))+(1:usfilter.dims.Mw(i)));
    end
    
    % calculate null spaces
    V = null([extplant.Cy extplant.Dyw]);
    W = null([extplant.Bu' extplant.Dzu']);
    V = blkdiag([V(1:extplant.dims.n-nxy,:) ; V(end-extplant.dims.mw+1:end,:)],eye(extplant.dims.mz));
    W = blkdiag([W(1:extplant.dims.n-nxu,:) ; W(end-extplant.dims.mz+1:end,:)],eye(extplant.dims.mw));
          
    % elaborate the actual synthesis LMIs
    if extplant.isdiscrete
        Zx = [[A11X;A21]'*X*[A11X;A21]-X11         [A11X;A21]'*X*extplant.Bw                Cz1'    ;
              ([A11X;A21]'*X*extplant.Bw)'         -Gamma_W+extplant.Bw'*X*extplant.Bw      Dstack' ;
              Cz1                                  Dstack                                   -Gamma_Z];

        Zy = [[A11Y A12]*Y*[A11Y A12]'-Y11         [A11Y A12]*Y*extplant.Cz'                Bw1     ;
              ([A11Y A12]*Y*extplant.Cz')'         -Gamma_Z+extplant.Cz*Y*extplant.Cz'      Dstack  ;
              Bw1'                                 Dstack'                                  -Gamma_W];
         
        PXY = [X eye(extplant.dims.n) ; eye(extplant.dims.n) Y];
    else
        Zx = [He(X11*A11X+X12*A21)        [X11 X12]*extplant.Bw     Cz1'    ;
              ([X11 X12]*extplant.Bw)'    -Gamma_W                  Dstack' ;
              Cz1                         Dstack                    -Gamma_Z];

        Zy = [He(A11Y*Y11+A12*Y21)        (extplant.Cz*[Y11;Y21])'     Bw1     ;
              (extplant.Cz*[Y11;Y21])     -Gamma_Z                     Dstack  ;
              Bw1'                        Dstack'                      -Gamma_W];
          
        IO = [eye(extplant.dims.n-nxy) zeros(extplant.dims.n-nxy,nxy)]';
        OI = [zeros(extplant.dims.n-nxu,nxu) eye(extplant.dims.n-nxu)];
        if nxu==0; PXY = [Y11 IO ; IO' X11]; else; PXY = [Y11 OI ; OI' X11]; end
    end

    ZX = V'*Zx*V;
    ZY = W'*Zy*W;
    constr = (PXY >= 0) + (ZX <= 0) + (ZY <= 0);
    
    % formulate the objective
    c = alpha(alpha~=0);
    obj = gamma*c(:);
    
    msgout('<strong>Solving synthesis problem...</strong>\n\n');
    
    % solve the problem
    options = mergestruct(gammasolver,sdpsettings());
    sol = optimize(constr, obj, options);
    if strfind(sol.info,'Infeasible'); error('The optimization failed: the constraints were found to be infeasible.'); end
    
    % check solution
    eigvals = [eig(double(PXY)) ; -eig(double(ZX)) ; -eig(double(ZY))];
    assert(all(eigvals >= 0), ['The SDP solver (''' gammasolver.solver ''') returned an infeasible solution for the synthesis problem. The most negative eigenvalue is ' num2str(min(eigvals)) '.']);
    
    % save auxiliary variables that are also required for the controller reconstruction
    aux.Zohat = Zohat; 
    aux.Zi = Zi;
    aux.Gammai = Gammai;
    aux.Gammao = Gammao;
    aux.Tu = Tu;
    aux.Ty = Ty;
    
    % extract the variables
    synresults.X11 = double(X11);
    synresults.X12 = double(X12);
    synresults.X22 = double(X22);
    synresults.Y11 = double(Y11);
    synresults.Y21 = double(Y21);
    synresults.Y22 = double(Y22);
    synresults.gamma = double(gamma);
    synresults.Gamma_W = double(Gamma_W);
    synresults.Gamma_Z = double(Gamma_Z);
    synresults.Dstack = double(Dstack);
    
end

function K = reconstruction_ct(stabplant,usfilter,extplant,synresults,aux,controllersolver)

    msgout('<strong>Reconstructing the stabilizing controller part...</strong>\n');
    
    % relax gamma values with 0.5% to ease the reconstruction
    if ~isfield(controllersolver,'relax'); controllersolver.relax = 0.005; end
    Gamma_W0 = synresults.Gamma_W + controllersolver.relax*diag(diag(synresults.Gamma_W));
    Gamma_Z0 = synresults.Gamma_Z + controllersolver.relax*diag(diag(synresults.Gamma_Z));
    
    % reconstruct the Lyapunov matrix 
    nxy = size(synresults.X12,2);
    nxu = size(synresults.Y21,1);
    
    if nxu > 0
        [U,S,V] = svd(eye(size(synresults.Y11))-[synresults.Y11 ; synresults.Y21]'*synresults.X11(:,1:extplant.dims.n-nxu));
        N = U*sqrt(S);
        M1 = V*sqrt(S);
        M2 = (-N\([synresults.Y11 ; synresults.Y21]'*synresults.X11(:,extplant.dims.n-nxu+1:end)))';

        IO = [eye(extplant.dims.n-nxu) zeros(extplant.dims.n-nxu,nxu)]';
        PI = [synresults.X11 IO ; M1' M2' zeros(extplant.dims.n-nxu)];
        Ptransf = [synresults.X11 IO ; IO' synresults.Y11];
        P = inv(PI'\Ptransf/PI);
    else
        [U,S,V] = svd(eye(size(synresults.X11))-[synresults.X11 synresults.X12]*synresults.Y11(:,1:extplant.dims.n-nxy));
        M = U*sqrt(S);
        N1 = V*sqrt(S);
        N2 = (-M\([synresults.X11 synresults.X12]*synresults.Y11(:,extplant.dims.n-nxy+1:end)))';

        IO = [eye(extplant.dims.n-nxy) zeros(extplant.dims.n-nxy,nxy)]';
        PI = [synresults.Y11 IO ; N1' N2' zeros(extplant.dims.n-nxy)];
        Ptransf = [synresults.Y11 IO ; IO' synresults.X11];
        P = PI'\Ptransf/PI;
    end

    % obtain the optimal closed loop as a function of the controller's state-space matrices
    A0 = zeros(2*extplant.dims.n-nxy-nxu); A0(1:extplant.dims.n,1:extplant.dims.n) = extplant.A;
    B0 = [extplant.Bw ; zeros(extplant.dims.n-nxy-nxu,extplant.dims.mw)];
    C0 = [extplant.Cz zeros(extplant.dims.mz,extplant.dims.n-nxy-nxu)];
    D0 = synresults.Dstack;

    Z = [A0'*P+P*A0   P*B0          C0'      ;
         B0'*P        -Gamma_W0     D0'      ; 
         C0           D0            -Gamma_Z0];  

    Buhat = [zeros(extplant.dims.n,extplant.dims.n-nxy-nxu) extplant.Bu ; eye(extplant.dims.n-nxy-nxu) zeros(extplant.dims.n-nxy-nxu, extplant.dims.mu)];
    Dzuhat = [zeros(extplant.dims.mz,extplant.dims.n-nxy-nxu) extplant.Dzu];
    Cyhat = [zeros(extplant.dims.n-nxy-nxu,extplant.dims.n) eye(extplant.dims.n-nxy-nxu) ; extplant.Cy zeros(extplant.dims.my,extplant.dims.n-nxy-nxu)];
    Dywhat = [zeros(extplant.dims.n-nxy-nxu,extplant.dims.mw) ; extplant.Dyw];

    U = [P*Buhat ; zeros(length(synresults.Gamma_W),size(Dzuhat,2)) ; Dzuhat]';
    V = [Cyhat Dywhat zeros(size(Cyhat,1),size(synresults.Gamma_Z,1))];

    % solve for the controller parameters
    switch controllersolver.solver
        case 'basiclmi'
            theta = basiclmi(Z,U,V,'Xmin,Shift');
        case 'lmilab'
            setlmis([]);
            theta = lmivar(2,[extplant.dims.n-nxy-nxu+extplant.dims.mu extplant.dims.n-nxy-nxu+extplant.dims.my]);
            receq = newlmi;
            lmiterm([receq 1 1 0],Z);
            lmiterm([receq 1 1 theta],U',V,'s')
            LMIs = getlmis;
            T = evalc('[~, xopt] = feasp(LMIs);');  % evalc to hide ugly output
            if contains(T,'INFEASIBLE'); error('The reconstruction problem appears to be infeasible. The synthesis may have returned ill-conditioned results.'); end
            theta = dec2mat(LMIs,xopt,theta);
        otherwise % through YALMIP
            theta = sdpvar(dims.n-nxy-nxu+dims.mu, dims.n-nxy-nxu+dims.my,'full');
            constr = [Z + U'*theta*V + V'*theta'*U <= 0];
            options = mergestruct(controllersolver,sdpsettings());
            optimize(constr,[],options);
            theta = double(theta);
            assert(all(eig(Z + U'*theta*V + V'*theta'*U) <= 0), ['The SDP solver (''' controllersolver.solver ''') returned an infeasible solution for the reconstruction problem.']);
    end
    
    Ak0 = theta(1:extplant.dims.n-nxy-nxu,1:extplant.dims.n-nxy-nxu);
    Bk0 = theta(1:extplant.dims.n-nxy-nxu,extplant.dims.n-nxy-nxu+1:end);
    Ck0 = theta(extplant.dims.n-nxy-nxu+1:end,1:extplant.dims.n-nxy-nxu);
    Dk0 = theta(extplant.dims.n-nxy-nxu+1:end,extplant.dims.n-nxy-nxu+1:end);
    
    % undo the swapping of control inputs and measured outputs that was
    % required for reducing the controller order
    Bk0 = Bk0*aux.Ty;
    Ck0 = aux.Tu*Ck0;
    Dk0 = aux.Tu*Dk0*aux.Ty;
    
    msgout('<strong>Recombining the controller with the unstable filter dynamics...</strong>\n');
    
    % split up the stabilizing part into the required partitions
    fictmy = size(Bk0,2)-stabplant.dims.my; % number of 'fictitious' measurements 
    Aa = Ak0;
    Ba1 = Bk0(:,1:fictmy);
    Ba2 = Bk0(:,fictmy+1:end);
    Ca1 = Ck0(stabplant.dims.mu+1:end,:);
    Ca2 = Ck0(1:stabplant.dims.mu,:);
    Da11 = Dk0(stabplant.dims.mu+1:end,1:fictmy);
    Da12 = Dk0(stabplant.dims.mu+1:end,fictmy+1:end);
    Da21 = Dk0(1:stabplant.dims.mu,1:fictmy);
    Da22 = Dk0(1:stabplant.dims.mu,fictmy+1:end);
    
    % include the dynamics of unstable weighting filters
    Ak0 = [usfilter.Ao-aux.Zohat*[Da11 ; Da21]      -aux.Zohat*[Ca1 ; Ca2]      -(aux.Gammao-aux.Zohat*[Da12 ; Da22])*aux.Zi;
           Ba1                                      Aa                          -Ba2*aux.Zi                                 ;
           Da11                                     Ca1                         usfilter.Ai-Da12*aux.Zi                     ];
    Bk0 = [aux.Gammao-aux.Zohat*[Da12 ; Da22];
           Ba2                               ;
           Da12                              ];
    Ck0 = [Da21 Ca2 (aux.Gammai-Da22*aux.Zi)];
    Dk0 = [Da22];
    
    msgout('<strong>Compensating for feedthrough matrix...</strong>\n\n');
    
    % compensate for feedthrough matrix Dyu
    if norm(stabplant.Dyu,1) > 0
        if norm(Dk0,1) > 0
            Myuk = eye(stabplant.dims.my)+stabplant.Dyu*Dk0; 
            Mkyu = eye(stabplant.dims.mu)+Dk0*stabplant.Dyu;
            if svds(Myuk,1,'smallest') < 1e-6
                error('The controller could not be reconstructed due to an algebraic loop in the generalized plant.');
            else
                Ck = Mkyu\Ck0;
                Dk = Mkyu\Dk0;
                Ak = Ak0-Bk0*stabplant.Dyu*Ck;
                Bk = Bk0/Myuk;
            end
        else
            Ak = Ak0-Bk0*stabplant.Dyu*Ck0;
            Bk = Bk0; Ck = Ck0; Dk = Dk0;
        end
    else
         Ak = Ak0; Bk = Bk0; Ck = Ck0; Dk = Dk0;
    end
    
    % return the controller
    K = ss(Ak,Bk,Ck,Dk);

end

function K = reconstruction_dt(stabplant,usfilter,extplant,synresults,aux,controllersolver)

    msgout('<strong>Reconstructing the stabilizing controller part...</strong>\n');
    
    % relax gamma values with 0.5% (default) to ease the reconstruction
    if ~isfield(controllersolver,'relax'); controllersolver.relax = 0.005; end
    Gamma_W0 = synresults.Gamma_W + controllersolver.relax*diag(diag(synresults.Gamma_W));
    Gamma_Z0 = synresults.Gamma_Z + controllersolver.relax*diag(diag(synresults.Gamma_Z));
    
    % reconstruct the Lyapunov matrix and find the closed loop
    nxy = size(synresults.X12,2);
    nxu = size(synresults.Y21,1);
    
     if nxu > 0
        Y = [synresults.Y11 synresults.Y21' ; synresults.Y21 synresults.Y22];
        Z = Y-inv(synresults.X11); Z11 = Z(1:end-nxu,1:end-nxu); Z12 = Z(1:end-nxu,end-nxu+1:end); Z22 = Z(end-nxu+1:end,end-nxu+1:end);
        Y22 = Y(end-nxu+1:end,end-nxu+1:end)-Z22+Z12'*pinv(Z11)*Z12; Y(end-nxu+1:end,end-nxu+1:end) = Y22;
        [U,S,V] = svd(eye(extplant.dims.n)-Y*synresults.X11);
        N = U*sqrt(S); N = N(:,1:end-nxu);
        M = V*sqrt(S); M = M(:,1:end-nxu);
        Pi22 = -M\synresults.X11*N;
        P = inv([Y N ; N' Pi22]);
     else
        X = [synresults.X11 synresults.X12 ; synresults.X12' synresults.X22];
        Z = X-inv(synresults.Y11); Z11 = Z(1:end-nxy,1:end-nxy); Z12 = Z(1:end-nxy,end-nxy+1:end); Z22 = Z(end-nxy+1:end,end-nxy+1:end);
        X22 = X(end-nxy+1:end,end-nxy+1:end)-Z22+Z12'*pinv(Z11)*Z12; X(end-nxy+1:end,end-nxy+1:end) = X22;
        [U,S,V] = svd(eye(extplant.dims.n)-X*synresults.Y11);
        N = V*sqrt(S); N = N(:,1:end-nxy);
        M = U*sqrt(S); M = M(:,1:end-nxy);
        P22 = -N\synresults.Y11*M;
        P = [X M ; M' P22];
     end
    
    % obtain the optimal closed loop as a function of the controller's state-space matrices
    A0 = zeros(2*extplant.dims.n-nxy-nxu); A0(1:extplant.dims.n,1:extplant.dims.n) = extplant.A;
    B0 = [extplant.Bw ; zeros(extplant.dims.n-nxy-nxu,extplant.dims.mw)];
    C0 = [extplant.Cz zeros(extplant.dims.mz,extplant.dims.n-nxy-nxu)];
    D0 = synresults.Dstack;

    Z = [-inv(P)                          A0                                B0                             zeros(size(P,1),size(D0,2));
         A0'                              -P                                zeros(size(P,1),size(B0,2))    C0'                        ; 
         B0'                              zeros(size(C0,1),size(P,2))       -Gamma_W0                      D0'                        ;
         zeros(size(C0,1),size(P,2))      C0                                D0                            -Gamma_Z0                   ];  

    Buhat = [zeros(extplant.dims.n,extplant.dims.n-nxy-nxu) extplant.Bu ; eye(extplant.dims.n-nxy-nxu) zeros(extplant.dims.n-nxy-nxu, extplant.dims.mu)];
    Dzuhat = [zeros(extplant.dims.mz,extplant.dims.n-nxy-nxu) extplant.Dzu];
    Cyhat = [zeros(extplant.dims.n-nxy-nxu,extplant.dims.n) eye(extplant.dims.n-nxy-nxu) ; extplant.Cy zeros(extplant.dims.my,extplant.dims.n-nxy-nxu)];
    Dywhat = [zeros(extplant.dims.n-nxy-nxu,extplant.dims.mw) ; extplant.Dyw];

    U = [Buhat' zeros(size(Buhat,2),size(Cyhat,2)+size(Dywhat,2)) Dzuhat'];
    V = [zeros(size(Cyhat,1),size(Buhat,1)) Cyhat Dywhat zeros(size(Cyhat,1),size(Dzuhat,1))];

    % solve for the controller parameters
    switch controllersolver.solver
        case 'basiclmi'
            theta = basiclmi(Z,U,V,'Xmin,Shift');
        case 'lmilab'
            setlmis([]);
            theta = lmivar(2,[extplant.dims.n-nxy-nxu+extplant.dims.mu extplant.dims.n-nxy-nxu+extplant.dims.my]);
            receq = newlmi;
            lmiterm([receq 1 1 0],Z);
            lmiterm([receq 1 1 theta],U',V,'s')
            LMIs = getlmis;
            T = evalc('[~, xopt] = feasp(LMIs);');  % evalc to hide ugly output
            if contains(T,'INFEASIBLE'); error('The reconstruction problem appears to be infeasible. The synthesis may have returned ill-conditioned results.'); end
            theta = dec2mat(LMIs,xopt,theta);
        otherwise % through YALMIP
            theta = sdpvar(dims.n-nxy-nxu+dims.mu, dims.n-nxy-nxu+dims.my,'full');
            constr = [Z + U'*theta*V + V'*theta'*U <= 0];
            options = mergestruct(controllersolver,sdpsettings());
            optimize(constr,[],options);
            theta = double(theta);
            assert(all(eig(Z + U'*theta*V + V'*theta'*U) <= 0), ['The SDP solver (''' controllersolver.solver ''') returned an infeasible solution for the reconstruction problem.']);
    end
    
    Ak0 = theta(1:extplant.dims.n-nxy-nxu,1:extplant.dims.n-nxy-nxu);
    Bk0 = theta(1:extplant.dims.n-nxy-nxu,extplant.dims.n-nxy-nxu+1:end);
    Ck0 = theta(extplant.dims.n-nxy-nxu+1:end,1:extplant.dims.n-nxy-nxu);
    Dk0 = theta(extplant.dims.n-nxy-nxu+1:end,extplant.dims.n-nxy-nxu+1:end);
    
    % undo the swapping of control inputs and measured outputs that was
    % required for reducing the controller order
    Bk0 = Bk0*aux.Ty;
    Ck0 = aux.Tu*Ck0;
    Dk0 = aux.Tu*Dk0*aux.Ty;
    
    msgout('<strong>Recombining the controller with the unstable filter dynamics...</strong>\n');
    
    % split up the stabilizing part into the required partitions
    fictmy = size(Bk0,2)-stabplant.dims.my; % number of 'fictitious' measurements 
    Aa = Ak0;
    Ba1 = Bk0(:,1:fictmy);
    Ba2 = Bk0(:,fictmy+1:end);
    Ca1 = Ck0(stabplant.dims.mu+1:end,:);
    Ca2 = Ck0(1:stabplant.dims.mu,:);
    Da11 = Dk0(stabplant.dims.mu+1:end,1:fictmy);
    Da12 = Dk0(stabplant.dims.mu+1:end,fictmy+1:end);
    Da21 = Dk0(1:stabplant.dims.mu,1:fictmy);
    Da22 = Dk0(1:stabplant.dims.mu,fictmy+1:end);
    
    % include the dynamics of unstable weighting filters
    Ak0 = [usfilter.Ao-aux.Zohat*[Da11 ; Da21]      -aux.Zohat*[Ca1 ; Ca2]      -(aux.Gammao-aux.Zohat*[Da12 ; Da22])*aux.Zi;
           Ba1                                      Aa                          -Ba2*aux.Zi                                 ;
           Da11                                     Ca1                         usfilter.Ai-Da12*aux.Zi                     ];
    Bk0 = [aux.Gammao-aux.Zohat*[Da12 ; Da22];
           Ba2                               ;
           Da12                              ];
    Ck0 = [Da21 Ca2 (aux.Gammai-Da22*aux.Zi)];
    Dk0 = [Da22];
    
    msgout('<strong>Compensating for feedthrough matrix...</strong>\n\n');
    
    % compensate for feedthrough matrix Dyu
    if norm(stabplant.Dyu,1) > 0
        if norm(Dk0,1) > 0
            Myuk = eye(stabplant.dims.my)+stabplant.Dyu*Dk0; 
            Mkyu = eye(stabplant.dims.mu)+Dk0*stabplant.Dyu;
            if svds(Myuk,1,'smallest') < 1e-6
                error('The controller could not be reconstructed due to an algebraic loop in the generalized plant.');
            else
                Ck = Mkyu\Ck0;
                Dk = Mkyu\Dk0;
                Ak = Ak0-Bk0*stabplant.Dyu*Ck;
                Bk = Bk0/Myuk;
            end
        else
            Ak = Ak0-Bk0*stabplant.Dyu*Ck0;
            Bk = Bk0; Ck = Ck0; Dk = Dk0;
        end
    else
         Ak = Ak0; Bk = Bk0; Ck = Ck0; Dk = Dk0;
    end
    
    % return the controller
    K = ss(Ak,Bk,Ck,Dk);

end

function sym = He(X)
    sym = X+X'; % notice: not divided by 2!
end

function [plant, nxy, nxu, Ty, Tu, T] = detectFCOR(plant)
    % Calculates a state-space realization that easily allows to reduce the
    % controller order while retaining convexity (if possible), i.e.:
    %
    %  / dx = Ax*x + Bw*w  + Bu*u
    % |  y  = Cy*x + Dyw*w
    %  \ z  = Cz*x + Dz*w  + Dzu*u
    %
    % with:  Cy = [Cy1 Cy2]     Dyw = [Dyw]
    %             [ 0   I ]           [ 0 ]
    % 
    %        or
    %
    %        Bu = [Bu1 0]     Dzu = [Dzu 0]
    %             [Bu2 I] 
    %
    % Apart from the new state-space realization, the transformation
    % matrices Ty and Tu that were used to obtain the new realization
    % are returned as well. nxy + nxu is the number of controller 
    % states that could be eliminated without losing optimality.
    
    % check which is the most advantageous option
    nxy = rank([plant.Cy plant.Dyw])-rank(plant.Dyw);
    nxu = rank([plant.Bu ; plant.Dzu])-rank(plant.Dzu);
    if nxu > nxy; nxy = 0; else; nxu = 0; end
    nxu = 0; nxy = 0;
    
    % transform to canonical form for continuous time
    if ~plant.isdiscrete
        if nxy==0
            % swap rows of Dzu
            [~,swap] = sort(sum(plant.Dzu,1)==0); 
            Tu = eye(length(swap)); Tu = Tu(swap,:); 
            Bu = plant.Bu*Tu;

            % look for an appropriate state transformation
            Bu(abs(Bu)<sqrt(eps)) = 0;
            [aux, j] = rref([fliplr(Bu) eye(plant.dims.n)]);
            Tx = fliplr(aux(:,plant.dims.mu+1:end)')';
            piv = j(j<=plant.dims.mu);
            T = zeros(plant.dims.mu);
            for i=1:nxu
                 T(piv(i),i) = 1;
                 T(i,piv(i)) = 1;
            end
            notsw = setdiff(1:plant.dims.mu,piv(1:nxu));
            T(sub2ind([plant.dims.mu plant.dims.mu],notsw,notsw)) = 1;
            T = blkdiag(T,eye(plant.dims.mu-length(T)));
            Tu = Tu*T;
            Ty = eye(plant.dims.my);
        else
            % swap columns of Dyw
            [~,swap] = sort(sum(plant.Dyw,2)==0);
            Ty = eye(length(swap)); Ty = Ty(swap,:);
            Cy = Ty*plant.Cy;

            % look for an approriate state transformation
            Cy(abs(Cy)<sqrt(eps)) = 0;
            [aux, j] = rref([fliplr(Cy') eye(plant.dims.n)]);
            Tx = inv(fliplr(aux(:,plant.dims.my+1:end)'));
            piv = j(j<=plant.dims.my);
            T = eye(plant.dims.my);
            for i=1:nxy
                 T(piv(i),i) = 1;
                 T(i,piv(i)) = 1;
            end
            notsw = setdiff(1:plant.dims.my,piv(1:nxy));
            T(sub2ind([plant.dims.my plant.dims.my],notsw,notsw)) = 1;
            T = blkdiag(T,eye(plant.dims.my-length(T)));
            Ty = T*Ty;
            Tu = eye(size(plant.dims.mu));
        end
        
        % obtain new full state-space matrices
        plant.A = Tx*plant.A/Tx;
        plant.Bu = Tx*plant.Bu*Tu;
        plant.Bw = Tx*plant.Bw;
        plant.Cy = Ty*plant.Cy/Tx;
        plant.Dyw = Ty*plant.Dyw;
        plant.Cz = plant.Cz/Tx;
        plant.Dzu = plant.Dzu*Tu; 
    end
    
    % rebalance
    [~,~,T] = balreal(ss(plant.A,plant.Bw,plant.Cz,plant.Dzw));
    T = blkdiag(T(1:end-nxu-nxy,1:end-nxu-nxy),T(end-nxu-nxy+1:end,end-nxu-nxy+1:end));
    plant.A = T*plant.A/T;
    plant.Bu = T*plant.Bu;
    plant.Bw = T*plant.Bw;
    plant.Cy = plant.Cy/T;
    plant.Cz = plant.Cz/T;
    
end

function msgout(msg,varargin)
    if verLessThan('matlab','7.13') || (usejava('jvm') && ~feature('ShowFigureWindows'))
        % remove HTML tags in messages if they're not supported (old MATLAB versions or from terminal)
        msg = regexprep(msg,'<(.*?)>','');
        fprintf(msg,varargin{:});
    else
        fprintf(msg,varargin{:});
    end
end