function [K, gamma] = projlem_d(prob,opts)
% PROJLEM_D Solves the LMI formulation of the standard problem in descriptor form
%
% [K,gamma] = PROJLEM(prob) returns a descriptor realization K of the
% controller that solves the SDP formulation of the standard problem,
% i.e. for the plant Pio (after elimination of impulsive and unstable 
% weighting filters)

% This file is part of hinfcd.
% Copyright (c) 2019, Laurens Jacobs, MECO Research Team @ KU Leuven. 
% 
% hinfcd is free software: you can redistribute it and/or modify it under 
% the terms of the GNU Lesser General Public License as published by the 
% Free Software Foundation, version 3.
% 
% hinfcd is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for
% more details.
% 
% You should have received a copy of the GNU Lesser General Public License
% along with hinfcd. If not, see <https://www.gnu.org/licenses/>.

    assert(prob.Ts==0,'Discrete-time problems not supported yet.');
    
    % 1. Standard synthesis
    switch lower(opts.synthesis.gammasolver)
        case 'yalmip'
            sdpvars = ct_synthesis_YALMIP(prob,opts);
        case 'cvx'
            sdpvars = ct_synthesis_CVX(prob,opts);
        case 'lmilab'
            sdpvars = ct_synthesis_LMILAB(prob,opts);
        otherwise
            prob.watchdog.error('LMI parser for synthesis problem not defined.'); 
    end
    
    % 2. Minimize convex envelope of Lyapunov matrix rank
    if opts.synthesis.minrank 
        switch lower(opts.synthesis.lyapunovsolver)
            case 'yalmip'
                sdpvars = ct_lyapunovrank_YALMIP(prob,sdpvars,opts);
            case 'cvx'
                sdpvars = ct_lyapunovrank_CVX(prob,sdpvars,opts);
            case 'lmilab'
                sdpvars = ct_lyapunovrank_LMILAB(prob,sdpvars,opts);
            otherwise
                prob.watchdog.error('LMI parser for synthesis problem not defined.'); 
        end
    end
    
    % 3. Reconstruct the Lyapunov matrix
    sdpvars = ct_lyapunov_reconstruction(sdpvars,opts);
    
    % 4. Controller reconstruction
    K = ct_controller_reconstruction(sdpvars,opts);
    
    % 5. Return gamma
    gamma = sqrt(sdpvars.sqgamma); 

end

%% CONTINUOUS TIME - SYNTHESIS

function sdpvars = ct_synthesis_YALMIP(prob,opts)
   if opts.synthesis.minrank
       plant = hinfcd.standard.projlem.genplant_d(prob);
   else
       plant = hinfcd.standard.projlem.genplant_d(prob,'bal');
   end
    
    % helpful index variables
    iscstr = [prob.specs.weight]==0;
    isobj = ~iscstr;
    nobj = sum(isobj);
    
    % make symmetry and equality constraints explicit
        % calculate explicit parametrization
        SYMX = hinfcd.util.explicitsym(plant.E);
        SYMYt = hinfcd.util.explicitsym(plant.E');
        EQW = hinfcd.util.explicitnull(plant.E,prob.nw());
        EQZt = hinfcd.util.explicitnull(plant.E',prob.nw());

    % define main variables
        % X, Y, W, Z
        X = sdpvar(size(SYMX,2),1,'full');
        Y = sdpvar(size(SYMYt,2),1,'full');
        W = sdpvar(size(EQW,2),1,'full');
        Z = sdpvar(size(EQZt,2),1,'full');
        
        % transform for ease of use
        X = reshape(SYMX*X,[plant.nx(),plant.nx()]);
        Y = reshape(SYMYt*Y,[plant.nx(),plant.nx()])';
        W = reshape(EQW*W,[plant.nx(),plant.nw()]);
        Z = reshape(EQZt*Z,[plant.nx(),plant.nw()])';
        
        % gamma
        sqgamma = sdpvar(nobj,1,'full');
        
    % construct matrices with slack variables
        Gammaw = sdpvar(plant.nw(),plant.nw(),'symmetric');
        Gammaz = sdpvar(plant.nz(),plant.nz(),'symmetric');
        Dzw = sdpvar(plant.nz(),plant.nw(),'full');
        
        % set constant parts & replace optimization variables
        j = 1;
        for i=1:length(prob.specs)
        	if prob.specs(i).weight==0
                Gammaz(prob.specs(i).out,prob.specs(i).out) = eye(length(prob.specs(i).out)); 
            else
                Gammaz(prob.specs(i).out,prob.specs(i).out) = sqgamma(j)*eye(length(prob.specs(i).out));
                j = j+1;
            end
            Gammaw(prob.specs(i).in,prob.specs(i).in) = eye(length(prob.specs(i).in)); 
            Dzw(prob.specs(i).out,prob.specs(i).in) = plant.Dzw(prob.specs(i).out,prob.specs(i).in);
        end
        
    % formulate LMIs
        % STABILITY 
        STAB = plant.NE()'*(0.5*[plant.E'*X+X'*plant.E   2*plant.E'           ;
                                 2*plant.E               plant.E*Y'+Y*plant.E'])*plant.NE();

        % FIRST PERFORMANCE LMI
        X1 = X(:,1:(plant.nx()-plant.smax));
        PERF1 = plant.NX()'*[plant.AX()'*X1+X1'*plant.AX()      X1'*plant.Bw()+plant.AX()'*W            plant.CzX()';
                             plant.Bw()'*X1+W'*plant.AX()       W'*plant.Bw()+plant.Bw()'*W-Gammaw      Dzw'        ;
                             plant.CzX()                        Dzw                                     -Gammaz     ]*plant.NX();
        
        % SECOND PERFORMANCE LMI
        Y1 = Y([1:(plant.nx()-plant.smax-plant.rmax) (plant.nx()-plant.smax+1):(plant.nx())],:);
        PERF2 = plant.NY()'*[Y1*plant.AY()'+plant.AY()*Y1'      Y1*plant.Cz'        plant.BwY()+plant.AY()*Z';
                             plant.Cz*Y1'                       -Gammaz             Dzw+plant.Cz*Z'          ;
                             Z*plant.AY()'+plant.BwY()'         Dzw'+Z*plant.Cz'    -Gammaw                  ]*plant.NY(); 
                
   % optimize
   obj = [prob.specs(isobj).weight]*sqgamma;
   cstr = [0.5*(STAB+STAB') + opts.synthesis.zerotol*eye(size(STAB)) >= 0 ;
           0.5*(PERF1+PERF1') - opts.synthesis.zerotol*eye(size(PERF1)) <= 0;
           0.5*(PERF2+PERF2') - opts.synthesis.zerotol*eye(size(PERF2)) <= 0];
   optimize(cstr,obj,opts.yalmip); 
   
   % put into the output structure
   sdpvars.X = double(X); 
   sdpvars.Y = double(Y); 
   sdpvars.W1 = double(W); 
   sdpvars.Z1 = double(Z); 
   sdpvars.Gammaw = double(Gammaw);
   sdpvars.Gammaz = double(Gammaz);
   sdpvars.D = double(Dzw);
   sdpvars.sqgamma = double(sqgamma); 
   sdpvars.plant = plant; 
   
   % make sure there are no NaNs
   sdpvars.X(isnan(sdpvars.X)) = 0;
   sdpvars.Y(isnan(sdpvars.Y)) = 0;
   sdpvars.W1(isnan(sdpvars.W1)) = 0;
   sdpvars.Z1(isnan(sdpvars.Z1)) = 0;
   sdpvars.Gammaw(isnan(sdpvars.Gammaw)) = 0;
   sdpvars.Gammaz(isnan(sdpvars.Gammaz)) = 0;
   sdpvars.D(isnan(sdpvars.D)) = 0;
   
end

function sdpvars = ct_synthesis_CVX(prob, opts)
    if opts.synthesis.minrank
        plant = hinfcd.standard.projlem.genplant_d(prob);
    else
        plant = hinfcd.standard.projlem.genplant_d(prob,'bal');
    end
   
    % helpful index variables
    iscstr = [prob.specs.weight]==0;
    isobj = ~iscstr;
    nobj = sum(isobj);
    
    % make symmetry and equality constraints explicit
        % calculate explicit parametrization
        SYMX = hinfcd.util.explicitsym(plant.E);
        SYMYt = hinfcd.util.explicitsym(plant.E');
        EQW = hinfcd.util.explicitnull(plant.E,prob.nw());
        EQZt = hinfcd.util.explicitnull(plant.E',prob.nw());
        
    % start CVX environment    
    cvx_begin SDP
        % X, Y, W, Z
        variable x(size(SYMX,2),1)
        variable y(size(SYMYt,2),1)
        variable w(size(EQW,2),1)
        variable z(size(EQZt,2),1)
        X = reshape(SYMX*x,[plant.nx(),plant.nx()]);
        Y = reshape(SYMYt*y,[plant.nx(),plant.nx()])';
        W = reshape(EQW*w,[plant.nx(),plant.nw()]);
        Z = reshape(EQZt*z,[plant.nx(),plant.nw()])';
        
        % gamma
        variable sqgamma(nobj,1)

        % construct structured matrix variables 
        j = 1;
        Gammaw = []; Gammaz = []; Dzw = [];
        for i=1:length(prob.specs)
            % Gammaw
            eval(['variable slackw' num2str(i) '(size(Gammaw,1),length(prob.specs(i).in))']);
            eval(['Gammaw = [Gammaw slackw' num2str(i) '; slackw' num2str(i) ''' eye(length(prob.specs(i).in))];'])
            
            % Gammaz
            eval(['variable slackz' num2str(i) '(size(Gammaz,1),length(prob.specs(i).out))']);
            if prob.specs(i).weight==0 
                eval(['Gammaz = [Gammaz slackz' num2str(i) '; slackz' num2str(i) ''' eye(length(prob.specs(i).out))];'])
            else
                eval(['Gammaz = [Gammaz slackz' num2str(i) '; slackz' num2str(i) ''' sqgamma(j)*eye(length(prob.specs(i).out))];'])
                j = j+1;
            end
            
            % Dzw
            eval(['variable slackdh' num2str(i) '(size(Dzw,1),length(prob.specs(i).in))']);
            eval(['variable slackdv' num2str(i) '(length(prob.specs(i).out),size(Dzw,2))']);
            eval(['Dzw = [Dzw slackdh' num2str(i) '; slackdv' num2str(i) ' plant.Dzw(prob.specs(i).out,prob.specs(i).in)];'])
        end
        
        % formulate LMIs
            % STABILITY 
            STAB = plant.NE()'*([plant.E'*X   plant.E'  ;
                                 plant.E      plant.E*Y'])*plant.NE();

            % FIRST PERFORMANCE LMI
            X1 = X(:,1:(plant.nx()-plant.smax));
            PERF1 = plant.NX()'*[plant.AX()'*X1+X1'*plant.AX()      X1'*plant.Bw()+plant.AX()'*W            plant.CzX()';
                                 plant.Bw()'*X1+W'*plant.AX()       W'*plant.Bw()+plant.Bw()'*W-Gammaw      Dzw'        ;
                                 plant.CzX()                        Dzw                                     -Gammaz     ]*plant.NX();

            % SECOND PERFORMANCE LMI
            Y1 = Y([1:(plant.nx()-plant.smax-plant.rmax) (plant.nx()-plant.smax+1):(plant.nx())],:);
            PERF2 = plant.NY()'*[Y1*plant.AY()'+plant.AY()*Y1'      Y1*plant.Cz'        plant.BwY()+plant.AY()*Z';
                                 plant.Cz*Y1'                       -Gammaz             Dzw+plant.Cz*Z'          ;
                                 Z*plant.AY()'+plant.BwY()'         Dzw'+Z*plant.Cz'    -Gammaw                  ]*plant.NY(); 
        
        % optimize
        minimize([prob.specs(isobj).weight]*sqgamma)
        subject to
            0.5*(STAB+STAB') + opts.synthesis.zerotol*eye(size(STAB)) >= 0
            0.5*(PERF1+PERF1') - opts.synthesis.zerotol*eye(size(PERF1)) <= 0
            0.5*(PERF2+PERF2') - opts.synthesis.zerotol*eye(size(PERF2)) <= 0
            
         eval(['cvx_solver ' opts.cvx.solver]);
         eval(['cvx_solver_settings(' strjoin(opts.cvx.solver_settings,',') ')']);
    cvx_end
    
   % put into the output structure
   sdpvars.X = X; 
   sdpvars.Y = Y; 
   sdpvars.W1 = W; 
   sdpvars.Z1 = Z; 
   sdpvars.Gammaw = Gammaw;
   sdpvars.Gammaz = Gammaz;
   sdpvars.D = Dzw;
   sdpvars.sqgamma = sqgamma;
   sdpvars.plant = plant; 
    
end

function sdpvars = ct_synthesis_LMILAB(prob,opts)
    if opts.synthesis.minrank
        plant = hinfcd.standard.projlem.genplant_d(prob);
    else
        plant = hinfcd.standard.projlem.genplant_d(prob,'bal');
    end
    
    % helpful index variables
    iscstr = [prob.specs.weight]==0;
    isobj = ~iscstr;
    nobj = sum(isobj);
    
    % make symmetry and equality constraints explicit
        % calculate explicit parametrization
        SYMX = hinfcd.util.explicitsym(plant.E);
        SYMYt = hinfcd.util.explicitsym(plant.E');
        EQW = hinfcd.util.explicitnull(plant.E,prob.nw());
        EQZt = hinfcd.util.explicitnull(plant.E',prob.nw());
    
        % helper functions
        function SYMXi = SYMXi(i,c)
            SYMXi = reshape(SYMX(:,i),[plant.nx() plant.nx()]); 
            if nargin>1
                SYMXi = SYMXi(:,c);
            end
        end
        function SYMYti = SYMYti(i,c)
            SYMYti = reshape(SYMYt(:,i),[plant.nx() plant.nx()]);
            if nargin>1
                SYMYti = SYMYti(:,c);
            end
        end
    
    % define main variables
    setlmis([]);
        % X, Y, W, Z
        nvar = 0; 
        x = zeros(size(SYMX,2),1); for i=1:length(x); x(i) = lmivar(1,[1 1]); end
        y = zeros(size(SYMYt,2),1); for i=1:length(y); y(i) = lmivar(1,[1 1]); end
        w = zeros(size(EQW,2),1);  for i=1:length(w); w(i) = lmivar(1,[1 1]); end
        z = zeros(size(EQZt,2),1); for i=1:length(z); z(i) = lmivar(1,[1 1]); end
        
        % gamma
        [sqgamma,nvar,ssqgamma] = lmivar(2,[nobj,1]);
    
    % construct matrices with slack variables
        % helper function
        function struc = offdiagonal(Nz,Nw,full)
            struc = zeros(sum(Nz),sum(Nw));
            Nw = [0 cumsum(Nw)]; Nz = [0 cumsum(Nz)]; 
            for ii=1:length(Nw)-1
                struc((Nz(ii+1)+1):end,(Nw(ii)+1):Nw(ii+1)) = 1;
                if nargin>=3 && full
                    struc(1:Nz(ii),(Nw(ii)+1):Nw(ii+1)) = 1;
                end
            end
        end

        % matrix structures
        iw = find(offdiagonal(cellfun(@length, {prob.specs.in}),cellfun(@length, {prob.specs.in}))); 
        iz = find(offdiagonal(cellfun(@length, {prob.specs.out}),cellfun(@length, {prob.specs.out}))); 
        id = find(offdiagonal(cellfun(@length, {prob.specs.out}),cellfun(@length, {prob.specs.in}),1)); 
        
        % create slack variables 
        sGammaw = zeros(prob.nw()); sGammaw(iw) = (1:length(iw))+nvar; sGammaw = sGammaw+sGammaw';
        sGammaz = zeros(prob.nz()); sGammaz(iz) = (1:length(iz))+nvar+length(iw); sGammaz = sGammaz+sGammaz';
        sD = zeros(prob.nz(),prob.nw()); sD(id) = (1:length(id))+nvar+length(iw)+length(iz); 
        
        % append optimization variables
        j = 1;
        Gammazcst = zeros(prob.nz()); 
        for i=1:length(prob.specs)
        	if prob.specs(i).weight==0
                Gammazcst(prob.specs(i).out,prob.specs(i).out) = eye(length(prob.specs(i).out)); 
            else
                sGammaz(prob.specs(i).out,prob.specs(i).out) = ssqgamma(j)*eye(length(prob.specs(i).out));
                j = j+1;
        	end
        end
        
        % transform to LMI variables
        Gammaw = lmivar(3,sGammaw);
        Gammaz = lmivar(3,sGammaz);
        D = lmivar(3,sD);
        
        % constant performance feedthrough matrices
        Dcst = zeros(size(plant.Dzw)); 
        for i=1:length(prob.specs)
            Dcst(prob.specs(i).out,prob.specs(i).in) = plant.Dzw(prob.specs(i).out,prob.specs(i).in);
        end
        
    % formulate LMIs
        % STABILITY 
        STAB = newlmi(); 
        for i=1:length(x); lmiterm([-STAB,1,1,x(i)],0.5*plant.E'*SYMXi(i),1,'s'); end
        for i=1:length(y); lmiterm([-STAB,2,2,y(i)],0.5*plant.E*SYMYti(i),1,'s'); end
        lmiterm([-STAB,2,1,0],plant.E);
        lmiterm([-STAB,0,0,0],plant.NE());
        lmiterm([-STAB,1,1,0],opts.synthesis.zerotol);
        lmiterm([-STAB,2,2,0],opts.synthesis.zerotol);

        % FIRST PERFORMANCE LMI
        PERF1 = newlmi();
        for i=1:length(x)
            lmiterm([PERF1,1,1,x(i)],plant.AX()'*SYMXi(i,1:(plant.nx-plant.smax)),1,'s'); 
            lmiterm([PERF1,2,1,x(i)],plant.Bw'*SYMXi(i,1:(plant.nx-plant.smax)),1);
        end
        for i=1:length(w)
            lmiterm([PERF1,2,1,w(i)],reshape(EQW(:,i),[plant.nx(), plant.nw()])'*plant.AX(),1); 
            lmiterm([PERF1,2,2,w(i)],reshape(EQW(:,i),[plant.nx(), plant.nw()])'*plant.Bw,1,'s'); 
        end
        lmiterm([PERF1,3,1,0],plant.CzX());
        if ~isempty(D); lmiterm([PERF1,3,2,D],1,1); end; lmiterm([PERF1,3,2,0],Dcst);
        if ~isempty(Gammaw); lmiterm([PERF1,2,2,Gammaw],-1,1); end; lmiterm([PERF1,2,2,0],-1);
        if ~isempty(Gammaz); lmiterm([PERF1,3,3,Gammaz],-1,1); end; lmiterm([PERF1,3,3,0],-Gammazcst); 
        lmiterm([PERF1,0,0,0],plant.NX());
        lmiterm([PERF1,1,1,0],-opts.synthesis.zerotol);
        lmiterm([PERF1,2,2,0],-opts.synthesis.zerotol);
        lmiterm([PERF1,3,3,0],-opts.synthesis.zerotol);

        % SECOND PERFORMANCE LMI
        PERF2 = newlmi();
        for i=1:length(y)
            lmiterm([PERF2,1,1,y(i)],plant.AY()*SYMYti(i,[1:(plant.nx()-plant.smax-plant.rmax) (plant.nx()-plant.smax+1):(plant.nx())]),1,'s'); 
            lmiterm([PERF2,2,1,y(i)],plant.Cz*SYMYti(i,[1:(plant.nx()-plant.smax-plant.rmax) (plant.nx()-plant.smax+1):(plant.nx())]),1);
        end
        for i=1:length(z)
           lmiterm([PERF2,1,3,z(i)],plant.AY()*reshape(EQZt(:,i),[plant.nx(), plant.nw()]),1); 
           lmiterm([PERF2,2,3,z(i)],plant.Cz*reshape(EQZt(:,i),[plant.nx(), plant.nw()]),1); 
        end
        lmiterm([PERF2,1,3,0],plant.BwY());
        if ~isempty(D); lmiterm([PERF2,2,3,D],1,1); end; lmiterm([PERF2,2,3,0],Dcst);
        if ~isempty(Gammaz); lmiterm([PERF2,2,2,Gammaz],-1,1); end; lmiterm([PERF2,2,2,0],-Gammazcst);
        if ~isempty(Gammaw); lmiterm([PERF2,3,3,Gammaw],-1,1); end; lmiterm([PERF2,3,3,0],-1); 
        lmiterm([PERF2,0,0,0],plant.NY());
        lmiterm([PERF2,1,1,0],-opts.synthesis.zerotol);
        lmiterm([PERF2,2,2,0],-opts.synthesis.zerotol);
        lmiterm([PERF2,3,3,0],-opts.synthesis.zerotol);
                
   % optimize
   problem = getlmis(); 
   c = zeros(decnbr(problem),1);
   if any(isobj)
        c(ssqgamma) = prob.specs(isobj).weight; 
   else
        c = [];
   end
   if isnumeric(opts.lmilab) && isvector(opts.lmilab)
        if isempty(c)
            [~,xopt] = feasp(problem,opts.lmilab);
        else
            [~,xopt] = mincx(problem,c,opts.lmilab);
        end
   else
        warning('Unknown options for LMILAB. Ignoring them. Check the mincx() documentation for details.'); 
        [~,xopt] = mincx(problem,c);
   end
 
   % get results
   for i=1:length(x); x(i) = dec2mat(problem,xopt,x(i)); end
   for i=1:length(y); y(i) = dec2mat(problem,xopt,y(i)); end
   for i=1:length(w); w(i) = dec2mat(problem,xopt,w(i)); end
   for i=1:length(z); z(i) = dec2mat(problem,xopt,z(i)); end
   
   % put into the output structure
   sdpvars.X = reshape(SYMX*x(:),[plant.nx() plant.nx()]);
   sdpvars.Y = reshape(SYMYt*y(:),[plant.nx() plant.nx()])';
   sdpvars.W1 = reshape(EQW*w(:),[plant.nx(), plant.nw()]);
   sdpvars.Z1 = reshape(EQZt*z(:),[plant.nx(), plant.nw()])';
   if ~isempty(Gammaw); sdpvars.Gammaw = dec2mat(problem,xopt,Gammaw)+eye(plant.nw()); else; sdpvars.Gammaw = eye(plant.nw()); end
   if ~isempty(Gammaz); sdpvars.Gammaz = dec2mat(problem,xopt,Gammaz)+Gammazcst; else; sdpvars.Gammaz = Gammazcst; end
   if ~isempty(D); sdpvars.D = dec2mat(problem,xopt,D)+Dcst; else; sdpvars.D = Dcst; end
   if ~isempty(sqgamma); sdpvars.sqgamma = dec2mat(problem,xopt,sqgamma); else; sdpvars.sqgamma = []; end
   sdpvars.plant = plant; 
  
end

function sdpvars = ct_lyapunovrank_YALMIP(prob, sdpvars, opts)
    plant = hinfcd.standard.projlem.genplant_d(prob,'bal'); 
    
    % make symmetry and equality constraints explicit
        % calculate explicit parametrization
        SYMX = hinfcd.util.explicitsym(plant.E);
        SYMYt = hinfcd.util.explicitsym(plant.E');
        EQW = hinfcd.util.explicitnull(plant.E,prob.nw());
        EQZt = hinfcd.util.explicitnull(plant.E',prob.nw());

    % define main variables
        % X, Y, W, Z
        X = sdpvar(size(SYMX,2),1,'full');
        Y = sdpvar(size(SYMYt,2),1,'full');
        W = sdpvar(size(EQW,2),1,'full');
        Z = sdpvar(size(EQZt,2),1,'full');
        
        % transform for ease of use
        X = reshape(SYMX*X,[plant.nx(),plant.nx()]);
        Y = reshape(SYMYt*Y,[plant.nx(),plant.nx()])';
        W = reshape(EQW*W,[plant.nx(),plant.nw()]);
        Z = reshape(EQZt*Z,[plant.nx(),plant.nw()])';
        
    % construct matrices with slack variables
        Gammaw = sdpvar(plant.nw(),plant.nw(),'symmetric');
        Gammaz = sdpvar(plant.nz(),plant.nz(),'symmetric');
        Dzw = sdpvar(plant.nz(),plant.nw(),'full');
        
        % set constant parts & replace optimization variables
        j = 1;
        for i=1:length(prob.specs)
        	if prob.specs(i).weight==0
                Gammaz(prob.specs(i).out,prob.specs(i).out) = eye(length(prob.specs(i).out)); 
            else
                Gammaz(prob.specs(i).out,prob.specs(i).out) = opts.synthesis.relaxation*sdpvars.sqgamma(j)*eye(length(prob.specs(i).out));
                j = j+1;
            end
            Gammaw(prob.specs(i).in,prob.specs(i).in) = eye(length(prob.specs(i).in)); 
            Dzw(prob.specs(i).out,prob.specs(i).in) = plant.Dzw(prob.specs(i).out,prob.specs(i).in);
        end
        
    % formulate LMIs
        % STABILITY 
        STAB = plant.NE()'*(0.5*[plant.E'*X+X'*plant.E   2*plant.E'           ;
                                 2*plant.E               plant.E*Y'+Y*plant.E'])*plant.NE();

        % FIRST PERFORMANCE LMI
        X1 = X(:,1:(plant.nx()-plant.smax));
        PERF1 = plant.NX()'*[plant.AX()'*X1+X1'*plant.AX()      X1'*plant.Bw()+plant.AX()'*W            plant.CzX()';
                             plant.Bw()'*X1+W'*plant.AX()       W'*plant.Bw()+plant.Bw()'*W-Gammaw      Dzw'        ;
                             plant.CzX()                        Dzw                                     -Gammaz     ]*plant.NX();
        
        % SECOND PERFORMANCE LMI
        Y1 = Y([1:(plant.nx()-plant.smax-plant.rmax) (plant.nx()-plant.smax+1):(plant.nx())],:);
        PERF2 = plant.NY()'*[Y1*plant.AY()'+plant.AY()*Y1'      Y1*plant.Cz'        plant.BwY()+plant.AY()*Z';
                             plant.Cz*Y1'                       -Gammaz             Dzw+plant.Cz*Z'          ;
                             Z*plant.AY()'+plant.BwY()'         Dzw'+Z*plant.Cz'    -Gammaw                  ]*plant.NY(); 
                
   % optimize
   obj = trace(plant.E'*X+plant.E*Y');
   cstr = [0.5*(STAB+STAB') + opts.synthesis.zerotol*eye(size(STAB)) >= 0 ;
           0.5*(PERF1+PERF1') - opts.synthesis.zerotol*eye(size(PERF1)) <= 0;
           0.5*(PERF2+PERF2') - opts.synthesis.zerotol*eye(size(PERF2)) <= 0];
   optimize(cstr,obj,opts.yalmip); 
   
   % calculate constraint tolerance
   sdpvars.epsilonX = max(eig(double(PERF1))); 
   sdpvars.epsilonY = max(eig(double(PERF2))); 
   
   % put into the output structure
   sdpvars.plant = plant; 
   sdpvars.X = double(X); 
   sdpvars.Y = double(Y); 
   sdpvars.W1 = double(W); 
   sdpvars.Z1 = double(Z); 
   sdpvars.Gammaw = double(Gammaw);
   sdpvars.Gammaz = double(Gammaz);
   sdpvars.D = double(Dzw);
   sdpvars.sqgamma = opts.synthesis.relaxation*sdpvars.sqgamma; 
      
   % make sure there are no NaNs
   sdpvars.X(isnan(sdpvars.X)) = 0;
   sdpvars.Y(isnan(sdpvars.Y)) = 0;
   sdpvars.W1(isnan(sdpvars.W1)) = 0;
   sdpvars.Z1(isnan(sdpvars.Z1)) = 0;
   sdpvars.Gammaw(isnan(sdpvars.Gammaw)) = 0;
   sdpvars.Gammaz(isnan(sdpvars.Gammaz)) = 0;
   sdpvars.D(isnan(sdpvars.D)) = 0;
   
end

function sdpvars = ct_lyapunovrank_CVX(prob, sdpvars, opts)
   plant = hinfcd.standard.projlem.genplant_d(prob,'bal');
  
    % make symmetry and equality constraints explicit
        % calculate explicit parametrization
        SYMX = hinfcd.util.explicitsym(plant.E);
        SYMYt = hinfcd.util.explicitsym(plant.E');
        EQW = hinfcd.util.explicitnull(plant.E,prob.nw());
        EQZt = hinfcd.util.explicitnull(plant.E',prob.nw());
        
    % start CVX environment    
    cvx_begin SDP
        % X, Y, W, Z
        variable x(size(SYMX,2),1)
        variable y(size(SYMYt,2),1)
        variable w(size(EQW,2),1)
        variable z(size(EQZt,2),1)
        X = reshape(SYMX*x,[plant.nx(),plant.nx()]);
        Y = reshape(SYMYt*y,[plant.nx(),plant.nx()])';
        W = reshape(EQW*w,[plant.nx(),plant.nw()]);
        Z = reshape(EQZt*z,[plant.nx(),plant.nw()])';

        % construct structured matrix variables 
        j = 1;
        Gammaw = []; Gammaz = []; Dzw = [];
        for i=1:length(prob.specs)
            % Gammaw
            eval(['variable slackw' num2str(i) '(size(Gammaw,1),length(prob.specs(i).in))']);
            eval(['Gammaw = [Gammaw slackw' num2str(i) '; slackw' num2str(i) ''' eye(length(prob.specs(i).in))];'])
            
            % Gammaz
            eval(['variable slackz' num2str(i) '(size(Gammaz,1),length(prob.specs(i).out))']);
            if prob.specs(i).weight==0 
                eval(['Gammaz = [Gammaz slackz' num2str(i) '; slackz' num2str(i) ''' eye(length(prob.specs(i).out))];'])
            else
                eval(['Gammaz = [Gammaz slackz' num2str(i) '; slackz' num2str(i) ''' opts.synthesis.relaxation*sdpvars.sqgamma(j)*eye(length(prob.specs(i).out))];'])
                j = j+1;
            end
            
            % Dzw
            eval(['variable slackdh' num2str(i) '(size(Dzw,1),length(prob.specs(i).in))']);
            eval(['variable slackdv' num2str(i) '(length(prob.specs(i).out),size(Dzw,2))']);
            eval(['Dzw = [Dzw slackdh' num2str(i) '; slackdv' num2str(i) ' plant.Dzw(prob.specs(i).out,prob.specs(i).in)];'])
        end
        
        % formulate LMIs
            % STABILITY 
            STAB = plant.NE()'*([plant.E'*X   plant.E'  ;
                                 plant.E      plant.E*Y'])*plant.NE();

            % FIRST PERFORMANCE LMI
            X1 = X(:,1:(plant.nx()-plant.smax));
            PERF1 = plant.NX()'*[plant.AX()'*X1+X1'*plant.AX()      X1'*plant.Bw()+plant.AX()'*W            plant.CzX()';
                                 plant.Bw()'*X1+W'*plant.AX()       W'*plant.Bw()+plant.Bw()'*W-Gammaw      Dzw'        ;
                                 plant.CzX()                        Dzw                                     -Gammaz     ]*plant.NX();

            % SECOND PERFORMANCE LMI
            Y1 = Y([1:(plant.nx()-plant.smax-plant.rmax) (plant.nx()-plant.smax+1):(plant.nx())],:);
            PERF2 = plant.NY()'*[Y1*plant.AY()'+plant.AY()*Y1'      Y1*plant.Cz'        plant.BwY()+plant.AY()*Z';
                                 plant.Cz*Y1'                       -Gammaz             Dzw+plant.Cz*Z'          ;
                                 Z*plant.AY()'+plant.BwY()'         Dzw'+Z*plant.Cz'    -Gammaw                  ]*plant.NY(); 
        
        % optimize
        minimize(trace(plant.E'*X+plant.E*Y'))
        subject to
            0.5*(STAB+STAB') + opts.synthesis.zerotol*eye(size(STAB)) >= 0
            0.5*(PERF1+PERF1') - opts.synthesis.zerotol*eye(size(PERF1)) <= 0
            0.5*(PERF2+PERF2') - opts.synthesis.zerotol*eye(size(PERF2)) <= 0
            
         eval(['cvx_solver ' opts.cvx.solver]);
         eval(['cvx_solver_settings(' strjoin(opts.cvx.solver_settings,',') ')']);
    cvx_end
    
   % calculate constraint tolerance
   sdpvars.epsilonX = max(eig(PERF1)); 
   sdpvars.epsilonY = max(eig(PERF2)); 
    
   % put into the output structure
   sdpvars.plant = plant; 
   sdpvars.X = X; 
   sdpvars.Y = Y; 
   sdpvars.W1 = W; 
   sdpvars.Z1 = Z; 
   sdpvars.Gammaw = Gammaw;
   sdpvars.Gammaz = Gammaz;
   sdpvars.D = Dzw;
   sdpvars.sqgamma = opts.synthesis.relaxation*sdpvars.sqgamma;
    
end

function sdpvars = ct_lyapunovrank_LMILAB(prob, sdpvars, opts)
    plant = hinfcd.standard.projlem.genplant_d(prob,'bal'); 
    
    % make symmetry and equality constraints explicit
        % calculate explicit parametrization
        SYMX = hinfcd.util.explicitsym(plant.E);
        SYMYt = hinfcd.util.explicitsym(plant.E');
        EQW = hinfcd.util.explicitnull(plant.E,prob.nw());
        EQZt = hinfcd.util.explicitnull(plant.E',prob.nw());
    
        % helper functions
        function SYMXi = SYMXi(i,c)
            SYMXi = reshape(SYMX(:,i),[plant.nx(), plant.nx()]); 
            if nargin>1
                SYMXi = SYMXi(:,c);
            end
        end
        function SYMYti = SYMYti(i,c)
            SYMYti = reshape(SYMYt(:,i),[plant.nx(), plant.nx()]);
            if nargin>1
                SYMYti = SYMYti(:,c);
            end
        end
    
    % define main variables
    setlmis([]);
        % X, Y, W, Z
        nvar = 0; 
        x = zeros(size(SYMX,2),1); for i=1:length(x); x(i) = lmivar(1,[1 1]); end
        y = zeros(size(SYMYt,2),1); for i=1:length(y); [y(i),nvar] = lmivar(1,[1 1]); end
        w = zeros(size(EQW,2),1);  for i=1:length(w); w(i) = lmivar(1,[1 1]); end
        z = zeros(size(EQZt,2),1); for i=1:length(z); [z(i),nvar] = lmivar(1,[1 1]); end
    
    % construct matrices with slack variables
        % helper function
        function struc = offdiagonal(Nz,Nw,full)
            struc = zeros(sum(Nz),sum(Nw));
            Nw = [0 cumsum(Nw)]; Nz = [0 cumsum(Nz)]; 
            for ii=1:length(Nw)-1
                struc((Nz(ii+1)+1):end,(Nw(ii)+1):Nw(ii+1)) = 1;
                if nargin>=3 && full
                    struc(1:Nz(ii),(Nw(ii)+1):Nw(ii+1)) = 1;
                end
            end
        end

        % matrix structures
        iw = find(offdiagonal(cellfun(@length, {prob.specs.in}),cellfun(@length, {prob.specs.in}))); 
        iz = find(offdiagonal(cellfun(@length, {prob.specs.out}),cellfun(@length, {prob.specs.out}))); 
        id = find(offdiagonal(cellfun(@length, {prob.specs.out}),cellfun(@length, {prob.specs.in}),1)); 
        
        % create slack variables 
        sGammaw = zeros(prob.nw()); sGammaw(iw) = (1:length(iw))+nvar; sGammaw = sGammaw+sGammaw';
        sGammaz = zeros(prob.nz()); sGammaz(iz) = (1:length(iz))+nvar+length(iw); sGammaz = sGammaz+sGammaz';
        sD = zeros(prob.nz(),prob.nw()); sD(id) = (1:length(id))+nvar+length(iw)+length(iz); 
        
        % append relaxed optimal objective values
        j = 1;
        Gammazcst = zeros(prob.nz()); 
        for i=1:length(prob.specs)
        	if prob.specs(i).weight==0
                Gammazcst(prob.specs(i).out,prob.specs(i).out) = eye(length(prob.specs(i).out)); 
            else
                Gammazcst(prob.specs(i).out,prob.specs(i).out) = opts.synthesis.relaxation*sdpvars.sqgamma(j)*eye(length(prob.specs(i).out));
                j = j+1;
        	end
        end
        
        % transform to LMI variables
        Gammaw = lmivar(3,sGammaw);
        Gammaz = lmivar(3,sGammaz);
        D = lmivar(3,sD);
        
        % constant performance feedthrough matrices
        Dcst = zeros(size(plant.Dzw)); 
        for i=1:length(prob.specs)
            Dcst(prob.specs(i).out,prob.specs(i).in) = plant.Dzw(prob.specs(i).out,prob.specs(i).in);
        end
        
    % formulate LMIs
        % STABILITY 
        STAB = newlmi(); 
        for i=1:length(x); lmiterm([-STAB,1,1,x(i)],0.5*plant.E'*SYMXi(i),1,'s'); end
        for i=1:length(y); lmiterm([-STAB,2,2,y(i)],0.5*plant.E*SYMYti(i),1,'s'); end
        lmiterm([-STAB,2,1,0],plant.E);
        lmiterm([-STAB,0,0,0],plant.NE());
        lmiterm([-STAB,1,1,0],opts.synthesis.zerotol);
        lmiterm([-STAB,2,2,0],opts.synthesis.zerotol);
        
        % FIRST PERFORMANCE LMI
        PERF1 = newlmi();
        for i=1:length(x)
            lmiterm([PERF1,1,1,x(i)],plant.AX()'*SYMXi(i,1:(plant.nx-plant.smax)),1,'s'); 
            lmiterm([PERF1,2,1,x(i)],plant.Bw'*SYMXi(i,1:(plant.nx-plant.smax)),1);
        end
        for i=1:length(w)
            lmiterm([PERF1,2,1,w(i)],reshape(EQW(:,i),[plant.nx(), plant.nw()])'*plant.AX(),1); 
            lmiterm([PERF1,2,2,w(i)],reshape(EQW(:,i),[plant.nx(), plant.nw()])'*plant.Bw,1,'s'); 
        end
        lmiterm([PERF1,3,1,0],plant.CzX());
        if ~isempty(D); lmiterm([PERF1,3,2,D],1,1); end; lmiterm([PERF1,3,2,0],Dcst);
        if ~isempty(Gammaw); lmiterm([PERF1,2,2,Gammaw],-1,1); end; lmiterm([PERF1,2,2,0],-1);
        if ~isempty(Gammaz); lmiterm([PERF1,3,3,Gammaz],-1,1); end; lmiterm([PERF1,3,3,0],-Gammazcst); 
        lmiterm([PERF1,0,0,0],plant.NX());
        lmiterm([PERF1,1,1,0],-opts.synthesis.zerotol);
        lmiterm([PERF1,2,2,0],-opts.synthesis.zerotol);
        lmiterm([PERF1,3,3,0],-opts.synthesis.zerotol);
        
        % SECOND PERFORMANCE LMI
        PERF2 = newlmi();
        for i=1:length(y)
            lmiterm([PERF2,1,1,y(i)],plant.AY()*SYMYti(i,[1:(plant.nx()-plant.smax-plant.rmax) (plant.nx()-plant.smax+1):(plant.nx())]),1,'s'); 
            lmiterm([PERF2,2,1,y(i)],plant.Cz*SYMYti(i,[1:(plant.nx()-plant.smax-plant.rmax) (plant.nx()-plant.smax+1):(plant.nx())]),1);
        end
        for i=1:length(z)
           lmiterm([PERF2,1,3,z(i)],plant.AY()*reshape(EQZt(:,i),[plant.nx(), plant.nw()]),1); 
           lmiterm([PERF2,2,3,z(i)],plant.Cz*reshape(EQZt(:,i),[plant.nx(), plant.nw()]),1); 
        end
        lmiterm([PERF2,1,3,0],plant.BwY());
        if ~isempty(D); lmiterm([PERF2,2,3,D],1,1); end; lmiterm([PERF2,2,3,0],Dcst);
        if ~isempty(Gammaz); lmiterm([PERF2,2,2,Gammaz],-1,1); end; lmiterm([PERF2,2,2,0],-Gammazcst);
        if ~isempty(Gammaw); lmiterm([PERF2,3,3,Gammaw],-1,1); end; lmiterm([PERF2,3,3,0],-1); 
        lmiterm([PERF2,0,0,0],plant.NY());
        lmiterm([PERF2,1,1,0],-opts.synthesis.zerotol);
        lmiterm([PERF2,2,2,0],-opts.synthesis.zerotol);
        lmiterm([PERF2,3,3,0],-opts.synthesis.zerotol);
                
    % minimize the trace (= convex envelope of rank) of stability LMI   
    problem = getlmis(); 
    c = zeros(decnbr(problem),1);
    for i=1:length(x); c(i) = trace(plant.E'*SYMXi(i)); end
    for i=1:length(y); c(i+length(x)) = trace(plant.E*SYMYti(i)); end
    
    if isnumeric(opts.lmilab) && isvector(opts.lmilab)
        [~,xopt] = mincx(problem,c,opts.lmilab);
    else
        warning('Unknown options for LMILAB. Ignoring them. Check the mincx() documentation for details.'); 
        [~,xopt] = mincx(problem,c);
    end

   % get results
   for i=1:length(x); x(i) = dec2mat(problem,xopt,x(i)); end
   for i=1:length(y); y(i) = dec2mat(problem,xopt,y(i)); end
   for i=1:length(w); w(i) = dec2mat(problem,xopt,w(i)); end
   for i=1:length(z); z(i) = dec2mat(problem,xopt,z(i)); end
   
   % calculate constraint tolerance
   PERF1 = showlmi(evallmi(problem,xopt),PERF1);
   PERF2 = showlmi(evallmi(problem,xopt),PERF2);
   sdpvars.epsilonX = max(eig(PERF1)); 
   sdpvars.epsilonY = max(eig(PERF2)); 
   
   % put into the output structure
   sdpvars.plant = plant;
   sdpvars.X = reshape(SYMX*x(:),[plant.nx(), plant.nx()]);
   sdpvars.Y = reshape(SYMYt*y(:),[plant.nx(), plant.nx()])';
   sdpvars.W1 = reshape(EQW*w(:),[plant.nx(), plant.nw()]);
   sdpvars.Z1 = reshape(EQZt*z(:),[plant.nx(), plant.nw()])';
   if ~isempty(Gammaw); sdpvars.Gammaw = dec2mat(problem,xopt,Gammaw)+eye(plant.nw()); else; sdpvars.Gammaw = eye(plant.nw()); end
   if ~isempty(Gammaz); sdpvars.Gammaz = dec2mat(problem,xopt,Gammaz)+Gammazcst; else; sdpvars.Gammaz = Gammazcst; end
   if ~isempty(D); sdpvars.D = dec2mat(problem,xopt,D)+Dcst; else; sdpvars.D = Dcst; end
   sdpvars.sqgamma = opts.synthesis.relaxation*sdpvars.sqgamma;
   
end

%% CONTINUOUS TIME - LYAPUNOV MATRIX RECONSTRUCTION 

function sdpvars = ct_lyapunov_reconstruction(sdpvars,opts)

    % get plant & synthesis data
    plant = sdpvars.plant; 
    n = plant.nx();
    re = rank(plant.E); 
    if opts.reconstruction.reducedorder
        nc1 = re-plant.r-plant.s; 
    else
        nc1 = re;
    end
    X = sdpvars.X;
    Y = sdpvars.Y; 
    
    % helper function
    function [M,N] = factorize(T,s0)
        [M,S,N] = svd(T); 
        if nargin<2; s0 = length(S); end
        ss = diag(S); s0 = min([s0,length(S)-sum(ss < opts.reconstruction.singvaltol)]); 
        M = M*sqrt(S(:,1:s0)); N = N*sqrt(S(:,1:s0)); 
    end

% plant is already in SVD canonical form!
%     % forward transformation (see Rehm & Allgower, ACC 2001)
%     [UE,SE,VE] = svd(plant.E); 
%     UE = UE*blkdiag(sqrt(SE(1:re,1:re)),eye(n-re)); 
%     VE = VE*blkdiag(sqrt(SE(1:re,1:re)),eye(n-re)); 
%     X = UE'*X/VE';
%     Y = VE'*Y/UE';
    
    % partition X and Y 
    X11 = X(1:re,1:re); 
    X21 = X((re+1):n,1:re);
    X22 = X((re+1):n,(re+1):n);
    Y11 = Y(1:re,1:re); 
    Y12 = Y(1:re,(re+1):n);
    Y22 = Y((re+1):n,(re+1):n);
    
    % iterative procedure to improve the low-rank approximation (see Xin, CDC 2003 & Xin et al., IEEE TAC 53(5) 2008)
    if nc1 < re
        if ~plant.i % the procedure described in the paper
            N1 = plant.NY()'*[eye(re) ; zeros(n-re+plant.nw()+plant.nz(),re)];
            J1 = plant.NY()'*[plant.A(:,1:re) ; plant.Cz(:,1:re) ; zeros(plant.nw(),re)];
            if isinf(plant.lambda); M1 = -N1; else; M1 = -plant.lambda*N1+J1; end
            [~,~,V] = svd(M1);
            iY11 = inv(Y11); 
            
            for i=1:opts.reconstruction.maxiter
                Z = V'*(X11-iY11)*V;
                Z11 = Z(1:nc1,1:nc1); 
                Z12 = Z(1:nc1,(nc1+1):re); 
                Z22 = Z((nc1+1):re,(nc1+1):re); 
                X11 = real(X11-V*blkdiag(zeros(nc1),Z22-Z12'*lsqminnorm(Z11,Z12))*V');
                s = svd(X11-iY11); 
                if max(s((nc1+1):re)) <= eps
                    break; 
                end
            end

            % update X in the SDP variables
            sdpvars.X = [X11 zeros(re,n-re) ; X21 X22]; 

        else % the dual procedure
            N1 = plant.NX()'*[eye(re) ; zeros(n-re+plant.nw()+plant.nz(),re)];
            J1 = plant.NX()'*[plant.A(1:re,:)' ; plant.Bw(1:re,:)' ; zeros(plant.nz(),re)];
            if isempty(plant.lambda); plant.lambda = Inf; end
            if isinf(plant.lambda); M1 = -N1; else; M1 = -plant.lambda*N1+J1; end
            [~,~,V] = svd(M1);
            iX11 = inv(X11); 
            
            for i=1:opts.reconstruction.maxiter
                Z = V'*(Y11-iX11)*V;
                Z11 = Z(1:nc1,1:nc1); 
                Z12 = Z(1:nc1,(nc1+1):re); 
                Z22 = Z((nc1+1):re,(nc1+1):re); 
                Y11 = real(Y11-V*blkdiag(zeros(nc1),Z22-Z12'*lsqminnorm(Z11,Z12))*V');
                s = svd(Y11-iX11); 
                if max(s((nc1+1):re)) <= eps
                    break; 
                end
            end

            % update Y in the SDP variables
            sdpvars.Y = [Y11 Y12 ; zeros(n-re,re) Y22]; 
        end
    end
    
    % obtain Lyapunov matrix 
        % M, U
        [M11,U11] = factorize(eye(re)-X11*Y11',nc1);
        [M22,U22] = factorize(eye(n-re)-X22*Y22');
        s = [kron(eye(re),M22) kron(U11,eye(n-re))]\vec(-X21*Y11'-X22*Y12');
        U12 = reshape(s(1:(re*(n-re))),[(n-re),re])';
        M21 = reshape(s((re*(n-re)+1):end),[n-re,nc1]);
        M = [M11 zeros(re,n-re) ; M21 M22];
        U = [U11 U12 ; zeros(n-re,nc1) U22];
        
        % N, R
        N11 = M11';
        R11 = -N11*Y11'/U11';
        R21 = zeros(n-re,nc1);
        R22 = eye(n-re);
        N22 = (eye(n-re)-U22')/Y22';
        N21 = (-N22*Y12'-R21*U11'-R22*U12')/Y11';
        N = [N11 zeros(nc1,n-re) ; N21 N22];
        R = [R11 zeros(nc1,n-re) ; R21 R22]; 

% plant is already in SVD canonical form!
%     % backward transformation
%     UEr = UE([1:nc1,(re+1):n],[1:nc1,(re+1):n]);
%     VEr = VE([1:nc1,(re+1):n],[1:nc1,(re+1):n]);
%     M = UE'\M*VEr';
%     N = UEr'\N*VE';
%     R = UEr'\R*VEr';
   
    % reconstruction of L, W2 and Z2
    sdpvars.L = diag([ones(nc1,1); zeros(n-re,1)]);
    sdpvars.Z2 = (-sdpvars.W1'-sdpvars.Z1*sdpvars.X')/M';
    sdpvars.W2 = -N*sdpvars.Z1'-R*sdpvars.Z2';
    sdpvars.P = [sdpvars.X M;
                 N         R];
             
    % useful for the reconstruction
    sdpvars.N = N;
    sdpvars.U = U; 
   
end

%% CONTINUOUS TIME - RECONSTRUCTION

function K = ct_controller_reconstruction(sdpvars,opts)

    plant = sdpvars.plant; 
     
    if strcmpi(opts.reconstruction.method,'transform')        
            warning(['Reconstruction of the transformed controller ' ...
                'parameters is not possible in the descriptor form ' ...
                'implementation. Directly reconstructing the actual ' ...
                'controller parameters instead.']); 
    end
    
    % obtain the required matrices
    P = sdpvars.P;
    W = [sdpvars.W1 ; sdpvars.W2];
    nc = length(P)-length(plant.E); 

    % reconstruct R
    Bue = [zeros(plant.nx(),nc)    plant.Bu;
           eye(nc)                 zeros(nc,plant.nu())]; 
    Dzue = [zeros(plant.nz(),nc)   plant.Dzu];
    R = [Bue'*P	Bue'*W Dzue'];

    % reconstruct S
    Cye = [zeros(nc,plant.nx())     eye(nc); 
           plant.Cy                 zeros(plant.ny(),nc)];
    Dywe = [zeros(nc,plant.nw()) ; plant.Dyw];
    S = [Cye Dywe zeros(nc+plant.ny(),plant.nz())];

    % reconstruct T
    A0 = blkdiag(plant.A,zeros(nc)); 
    B0 = [plant.Bw ; zeros(nc,plant.nw())];
    C0 = [plant.Cz zeros(plant.nz(),nc)];
    D0 = sdpvars.D;
    T = [A0'*P+P'*A0    P'*B0+A0'*W                 C0'            ;
         B0'*P+W'*A0    W'*B0+B0'*W-sdpvars.Gammaw  D0'            ;
         C0             D0                          -sdpvars.Gammaz];

    % solve
    switch lower(opts.reconstruction.solver)
        case 'lmilab'
            Theta = ct_reconstruction_LMILAB(R,S,T,opts);
        case 'yalmip'
            Theta = ct_reconstruction_YALMIP(R,S,T,opts); 
        case 'cvx'
            Theta = ct_reconstruction_CVX(R,S,T,opts); 
        case 'basiclmi'
            Theta = ct_reconstruction_BASICLMI(R,S,T); 
        otherwise
            error('LMI solver for feasibility problem (controller reconstruction) not found. Supported interfaces are ''lmilab'', ''yalmip'', ''cvx'' and ''basiclmi''.'); 
    end
    
    % transform back to controller for the original plant
    K = hinfcd.dss(Theta(1:nc,1:nc), ...
                   Theta(1:nc,(nc+1):end), ...
                   Theta((nc+1):end,1:nc), ...
                   Theta((nc+1):end,(nc+1):end), ...
                   sdpvars.L, plant.Ts);

end

function Theta = ct_reconstruction_YALMIP(R,S,T,opts)
    % solve for Theta
    Theta = sdpvar(size(R,1),size(S,1),'full'); 
    optimize([T+R'*Theta*S+S'*Theta'*R <= 0], [], opts.yalmip);
    Theta = double(Theta); 
end

function Theta = ct_reconstruction_CVX(R,S,T,opts)
    % solve for Theta
    cvx_begin SDP
        variable Theta(size(R,1), size(S,1))
        subject to
            T+R'*Theta*S+S'*Theta'*R <= 0
        eval(['cvx_solver ' opts.cvx.solver]);
        eval(['cvx_solver_settings(' strjoin(opts.cvx.solver_settings,',') ')']);
    cvx_end 
end

function Theta = ct_reconstruction_LMILAB(R,S,T,opts)
    % solve for Theta
    setlmis([]);
    Theta = lmivar(2,[size(R,1) size(S,1)]);
    Z = newlmi(); 
    lmiterm([Z,1,1,0],0.5*(T+T'));
    lmiterm([Z,1,1,Theta],R',S,'s');
    problem = getlmis();
    if isnumeric(opts.lmilab) && isvector(opts.lmilab)
        [~,xfeas] = feasp(problem,opts.lmilab,opts.reconstruction.zerotol);
    else
        warning('Unknown options for LMILAB. Ignoring them. Check the mincx() documentation for details.'); 
        [~,xfeas] = feasp(problem,zeros(1,5),opts.synthesis.zerotol);
    end
    Theta = dec2mat(problem,xfeas,Theta);
end

function Theta = ct_reconstruction_BASICLMI(R,S,T)
    Theta = basiclmi(T,R,S,'Xmin,Shift');
end
