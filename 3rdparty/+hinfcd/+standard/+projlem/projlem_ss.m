function [K, gamma] = projlem_ss(prob,opts)
% PROJLEM_SS Solves the LMI formulation of the standard problem in state-space form
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
    plant = hinfcd.standard.projlem.genplant_ss(prob);
    
    % helpful index variables
    iscstr = [prob.specs.weight]==0;
    isobj = ~iscstr;
    nobj = sum(isobj);
    
    % define main variables
        % X, Y
        X = sdpvar(plant.nx,plant.nx,'symmetric');
        Y = sdpvar(plant.nx,plant.nx,'symmetric');
        
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
        STAB = [X                eye(plant.nx()); 
                eye(plant.nx())  Y'             ];

        % FIRST PERFORMANCE LMI
        PERF1 = plant.NX()'*[plant.A'*X+X'*plant.A      X'*plant.Bw()   plant.Cz';
                             plant.Bw'*X                -Gammaw         Dzw'     ;
                             plant.Cz                   Dzw              -Gammaz ]*plant.NX();
        
        % SECOND PERFORMANCE LMI
        PERF2 = plant.NY()'*[Y*plant.A'+plant.A*Y'      Y*plant.Cz'     plant.Bw;
                             plant.Cz*Y'                -Gammaz         Dzw     ;
                             plant.Bw'                  Dzw'            -Gammaw ]*plant.NY(); 
                
   % optimize
   obj = vertcat(prob.specs(isobj).weight)'*sqgamma;
   cstr = [0.5*(STAB+STAB') + opts.synthesis.zerotol*eye(size(STAB)) >= 0 ;
           0.5*(PERF1+PERF1') - opts.synthesis.zerotol*eye(size(PERF1)) <= 0;
           0.5*(PERF2+PERF2') - opts.synthesis.zerotol*eye(size(PERF2)) <= 0];
   optimize(cstr,obj,opts.yalmip); 
   
   % put into the output structure
   sdpvars.X = double(X); 
   sdpvars.Y = double(Y); 
   sdpvars.Gammaw = double(Gammaw);
   sdpvars.Gammaz = double(Gammaz);
   sdpvars.D = double(Dzw);
   sdpvars.sqgamma = double(sqgamma); 
   sdpvars.plant = plant; 
   
end

function sdpvars = ct_synthesis_CVX(prob, opts)
    plant = hinfcd.standard.projlem.genplant_ss(prob);
   
    % helpful index variables
    iscstr = [prob.specs.weight]==0;
    isobj = ~iscstr;
    nobj = sum(isobj);
        
    % start CVX environment    
    cvx_begin SDP
        % X, Y
        variable X(size(plant.A)) symmetric
        variable Y(size(plant.A)) symmetric
        
        % gamma
        variable sqgamma(nobj,1)

        % construct structured matrix variables 
        j = 1;
        Gammaw = eye(length(prob.specs(1).in));
        Gammaz = sqgamma(j)*eye(length(prob.specs(1).out)); j = j+1;
        Dzw = plant.Dzw(prob.specs(1).out,prob.specs(1).in);
        for i=2:length(prob.specs)
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
            STAB = [X                eye(plant.nx()); 
                    eye(plant.nx())  Y'             ];

            % FIRST PERFORMANCE LMI
            PERF1 = plant.NX()'*[plant.A'*X+X'*plant.A      X'*plant.Bw     plant.Cz';
                                 plant.Bw'*X                -Gammaw         Dzw'     ;
                                 plant.Cz                   Dzw              -Gammaz ]*plant.NX();

            % SECOND PERFORMANCE LMI
            PERF2 = plant.NY()'*[Y*plant.A'+plant.A*Y'      Y*plant.Cz'     plant.Bw;
                                 plant.Cz*Y'                -Gammaz         Dzw     ;
                                 plant.Bw'                  Dzw'            -Gammaw ]*plant.NY(); 
                
        % optimize
        minimize([prob.specs(isobj).weight]*sqgamma)
        subject to
            0.5*(STAB+STAB') + opts.synthesis.zerotol*eye(size(STAB)) >= 0
            0.5*(PERF1+PERF1') - opts.synthesis.zerotol*eye(size(PERF1)) <= 0
            0.5*(PERF2+PERF2') - opts.synthesis.zerotol*eye(size(PERF2)) <= 0
            
         eval(['cvx_solver ' opts.cvx.solver]);
         eval(['cvx_precision ' opts.cvx.precision]);
         eval(['cvx_solver_settings(' strjoin(opts.cvx.solver_settings,',') ')']);
    cvx_end
    
   % put into the output structure
   sdpvars.X = X; 
   sdpvars.Y = Y; 
   sdpvars.Gammaw = Gammaw;
   sdpvars.Gammaz = Gammaz;
   sdpvars.D = Dzw;
   sdpvars.sqgamma = sqgamma;
   sdpvars.plant = plant; 
    
end

function sdpvars = ct_synthesis_LMILAB(prob,opts)
    plant = hinfcd.standard.projlem.genplant_ss(prob);
    
    % helpful index variables
    iscstr = [prob.specs.weight]==0;
    isobj = ~iscstr;
    nobj = sum(isobj);
    
    % define main variables
    setlmis([]);
        % X, Y
        X = lmivar(1,[length(plant.A),1]);
        Y = lmivar(1,[length(plant.A),1]); 
        
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
        lmiterm([-STAB,1,1,X],1,1); 
        lmiterm([-STAB,2,2,Y],1,1); 
        lmiterm([-STAB,2,1,0],1);
        lmiterm([-STAB,1,1,0],opts.synthesis.zerotol);
        lmiterm([-STAB,2,2,0],opts.synthesis.zerotol);

        % FIRST PERFORMANCE LMI
        PERF1 = newlmi();
        lmiterm([PERF1,1,1,X],plant.A',1,'s'); 
        lmiterm([PERF1,2,1,X],plant.Bw',1);
        lmiterm([PERF1,3,1,0],plant.Cz);
        if ~isempty(D); lmiterm([PERF1,3,2,D],1,1); end; lmiterm([PERF1,3,2,0],Dcst);
        if ~isempty(Gammaw); lmiterm([PERF1,2,2,Gammaw],-1,1); end; lmiterm([PERF1,2,2,0],-1);
        if ~isempty(Gammaz); lmiterm([PERF1,3,3,Gammaz],-1,1); end; lmiterm([PERF1,3,3,0],-Gammazcst); 
        lmiterm([PERF1,0,0,0],plant.NX());
        lmiterm([PERF1,1,1,0],-opts.synthesis.zerotol);
        lmiterm([PERF1,2,2,0],-opts.synthesis.zerotol);
        lmiterm([PERF1,3,3,0],-opts.synthesis.zerotol);

        % SECOND PERFORMANCE LMI
        PERF2 = newlmi();
        lmiterm([PERF2,1,1,Y],1,plant.A','s'); 
        lmiterm([PERF2,1,2,Y],1,plant.Cz',1);
        lmiterm([PERF2,1,3,0],plant.Bw);
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
 
   % put into the output structure
   sdpvars.X = dec2mat(problem,xopt,X);
   sdpvars.Y = dec2mat(problem,xopt,Y);
   if ~isempty(Gammaw); sdpvars.Gammaw = dec2mat(problem,xopt,Gammaw)+eye(plant.nw()); else; sdpvars.Gammaw = eye(plant.nw()); end
   if ~isempty(Gammaz); sdpvars.Gammaz = dec2mat(problem,xopt,Gammaz)+Gammazcst; else; sdpvars.Gammaz = Gammazcst; end
   if ~isempty(D); sdpvars.D = dec2mat(problem,xopt,D)+Dcst; else; sdpvars.D = Dcst; end
   if any(isobj); sdpvars.sqgamma = dec2mat(problem,xopt,sqgamma); else; sdpvars.sqgamma = []; end
   sdpvars.plant = plant; 
  
end

function sdpvars = ct_lyapunovrank_YALMIP(prob, sdpvars, opts)
    plant = sdpvars.plant; 
    
    % helpful index variables
    iscstr = [prob.specs.weight]==0;
    isobj = ~iscstr;
    nobj = sum(isobj);
    
    % define main variables
        % X, Y
        X = sdpvar(plant.nx,plant.nx,'symmetric');
        Y = sdpvar(plant.nx,plant.nx,'symmetric');
        
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
        STAB = [X                eye(plant.nx()); 
                eye(plant.nx())  Y'             ];

        % FIRST PERFORMANCE LMI
        PERF1 = plant.NX()'*[plant.A'*X+X'*plant.A      X'*plant.Bw()   plant.Cz';
                             plant.Bw'*X                -Gammaw         Dzw'     ;
                             plant.Cz                   Dzw              -Gammaz ]*plant.NX();
        
        % SECOND PERFORMANCE LMI
        PERF2 = plant.NY()'*[Y*plant.A'+plant.A*Y'      Y*plant.Cz'     plant.Bw;
                             plant.Cz*Y'                -Gammaz         Dzw     ;
                             plant.Bw'                  Dzw'            -Gammaw ]*plant.NY(); 
                
   % optimize
   obj = trace(X)+trace(Y');
   cstr = [0.5*(STAB+STAB') + opts.synthesis.zerotol*eye(size(STAB)) >= 0 ;
           0.5*(PERF1+PERF1') - opts.synthesis.zerotol*eye(size(PERF1)) <= 0;
           0.5*(PERF2+PERF2') - opts.synthesis.zerotol*eye(size(PERF2)) <= 0];
   optimize(cstr,obj,opts.yalmip); 
   
   % put into the output structure
   sdpvars.plant = plant; 
   sdpvars.X = double(X); 
   sdpvars.Y = double(Y); 
   sdpvars.Gammaw = double(Gammaw);
   sdpvars.Gammaz = double(Gammaz);
   sdpvars.D = double(Dzw);
   sdpvars.sqgamma = opts.synthesis.relaxation*sdpvars.sqgamma; 
   
end

function sdpvars = ct_lyapunovrank_CVX(prob, sdpvars, opts)
    plant = sdpvars.plant; 
        
    % start CVX environment    
    cvx_begin SDP
        % X, Y
        variable X(size(plant.A)) symmetric
        variable Y(size(plant.A)) symmetric

        % construct structured matrix variables 
        j = 1;
        Gammaw = eye(length(prob.specs(1).in));
        Gammaz = opts.synthesis.relaxation*sdpvars.sqgamma(j)*eye(length(prob.specs(1).out)); j = j+1;
        Dzw = plant.Dzw(prob.specs(1).out,prob.specs(1).in);
        for i=2:length(prob.specs)
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
            STAB = [X                eye(plant.nx()); 
                    eye(plant.nx())  Y'             ];

            % FIRST PERFORMANCE LMI
            PERF1 = plant.NX()'*[plant.A'*X+X'*plant.A      X'*plant.Bw     plant.Cz';
                                 plant.Bw'*X                -Gammaw         Dzw'     ;
                                 plant.Cz                   Dzw              -Gammaz ]*plant.NX();

            % SECOND PERFORMANCE LMI
            PERF2 = plant.NY()'*[Y*plant.A'+plant.A*Y'      Y*plant.Cz'     plant.Bw;
                                 plant.Cz*Y'                -Gammaz         Dzw     ;
                                 plant.Bw'                  Dzw'            -Gammaw ]*plant.NY(); 
                
        % optimize
        minimize(trace(X)+trace(Y'))
        subject to
            0.5*(STAB+STAB') + opts.synthesis.zerotol*eye(size(STAB)) >= 0
            0.5*(PERF1+PERF1') - opts.synthesis.zerotol*eye(size(PERF1)) <= 0
            0.5*(PERF2+PERF2') - opts.synthesis.zerotol*eye(size(PERF2)) <= 0
            
         eval(['cvx_solver ' opts.cvx.solver]);
         eval(['cvx_precision ' opts.cvx.precision]);
         eval(['cvx_solver_settings(' strjoin(opts.cvx.solver_settings,',') ')']);
    cvx_end
    
   % put into the output structure
   sdpvars.plant = plant; 
   sdpvars.X = X; 
   sdpvars.Y = Y; 
   sdpvars.Gammaw = Gammaw;
   sdpvars.Gammaz = Gammaz;
   sdpvars.D = Dzw;
   sdpvars.sqgamma = opts.synthesis.relaxation*sdpvars.sqgamma;
    
end

function sdpvars = ct_lyapunovrank_LMILAB(prob, sdpvars, opts)
    plant = sdpvars.plant; 
    
    % define main variables
    setlmis([]);
        nvar = 0; 
        % X, Y
        [X,nvar,sX] = lmivar(1,[length(plant.A),1]);
        [Y,nvar,sY] = lmivar(1,[length(plant.A),1]); 
    
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
        lmiterm([-STAB,1,1,X],1,1); 
        lmiterm([-STAB,2,2,Y],1,1); 
        lmiterm([-STAB,2,1,0],1);
        lmiterm([-STAB,1,1,0],opts.synthesis.zerotol);
        lmiterm([-STAB,2,2,0],opts.synthesis.zerotol);

        % FIRST PERFORMANCE LMI
        PERF1 = newlmi();
        lmiterm([PERF1,1,1,X],plant.A',1,'s'); 
        lmiterm([PERF1,2,1,X],plant.Bw',1);
        lmiterm([PERF1,3,1,0],plant.Cz);
        if ~isempty(D); lmiterm([PERF1,3,2,D],1,1); end; lmiterm([PERF1,3,2,0],Dcst);
        if ~isempty(Gammaw); lmiterm([PERF1,2,2,Gammaw],-1,1); end; lmiterm([PERF1,2,2,0],-1);
        if ~isempty(Gammaz); lmiterm([PERF1,3,3,Gammaz],-1,1); end; lmiterm([PERF1,3,3,0],-Gammazcst); 
        lmiterm([PERF1,0,0,0],plant.NX());
        lmiterm([PERF1,1,1,0],-opts.synthesis.zerotol);
        lmiterm([PERF1,2,2,0],-opts.synthesis.zerotol);
        lmiterm([PERF1,3,3,0],-opts.synthesis.zerotol);
        
        % SECOND PERFORMANCE LMI
        PERF2 = newlmi(); 
        lmiterm([PERF2,1,1,Y],1,plant.A','s'); 
        lmiterm([PERF2,1,2,Y],1,plant.Cz',1);
        lmiterm([PERF2,1,3,0],plant.Bw);
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
    c([diag(sX);diag(sY)]) = ones(2*length(plant.A),1);
    
    if isnumeric(opts.lmilab) && isvector(opts.lmilab)
        [~,xopt] = mincx(problem,c,opts.lmilab);
    else
        warning('Unknown options for LMILAB. Ignoring them. Check the mincx() documentation for details.'); 
        [~,xopt] = mincx(problem,c);
    end

   % put into the output structure
   sdpvars.plant = plant;
   sdpvars.X = dec2mat(problem,xopt,X);
   sdpvars.Y = dec2mat(problem,xopt,Y);
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
    X = sdpvars.X;
    Y = sdpvars.Y; 
    s = svd(eye(size(X))-X*Y);
    stol = length(s)-sum(s < opts.reconstruction.singvaltol);
    if opts.reconstruction.reducedorder
        nc1 = min(n-plant.r-plant.s,stol); 
    else
        nc1 = n;
    end
%    nc1 = n-1;
    
    % helper function
    function [M,N] = factorize(X,s0)
        [M,S,N] = svd(X); 
        if nargin<2; s0 = length(S); end
        M = M*sqrt(S(:,1:s0)); N = N*sqrt(S(:,1:s0)); 
    end

    % iterative procedure to improve the low-rank approximation (see Xin, CDC 2003 & Xin et al., IEEE TAC 53(5) 2008)
    if nc1 < n
        if ~plant.i % the procedure described in the paper
            N1 = plant.NY()'*[eye(n) ; zeros(plant.nw()+plant.nz(),n)];
            J1 = plant.NY()'*[plant.A(:,1:n) ; plant.Cz(:,1:n) ; zeros(plant.nw(),n)];
            if isinf(plant.lambda); M1 = -N1; else; M1 = -plant.lambda*N1+J1; end
            [~,~,V] = svd(M1);
            iY = inv(Y); 
            
            for i=1:opts.reconstruction.maxiter
                Z = V'*(X-iY)*V;
                Z11 = Z(1:nc1,1:nc1); 
                Z12 = Z(1:nc1,(nc1+1):n); 
                Z22 = Z((nc1+1):n,(nc1+1):n); 
                X = real(X-V*blkdiag(zeros(nc1),Z22-Z12'*lsqminnorm(Z11,Z12))*V');
                s = svd(X-iY); 
                if max(s((nc1+1):n)) <= eps
                    break; 
                end
            end

            % update X in the SDP variables
            sdpvars.X = X;

        else % the dual procedure
            N1 = plant.NX()'*[eye(n) ; zeros(plant.nw()+plant.nz(),n)];
            J1 = plant.NX()'*[plant.A(1:n,:)' ; plant.Bw(1:n,:)' ; zeros(plant.nz(),n)];
            if isempty(plant.lambda); plant.lambda = Inf; end
            if isinf(plant.lambda); M1 = -N1; else; M1 = -plant.lambda*N1+J1; end
            [~,~,V] = svd(M1);
            iX = inv(X); 
            
            for i=1:opts.reconstruction.maxiter
                Z = V'*(Y-iX)*V;
                Z11 = Z(1:nc1,1:nc1); 
                Z12 = Z(1:nc1,(nc1+1):n); 
                Z22 = Z((nc1+1):n,(nc1+1):n); 
                Y = real(Y-V*blkdiag(zeros(nc1),Z22-Z12'*lsqminnorm(Z11,Z12))*V');
                s = svd(Y-iX); 
                if max(s((nc1+1):n)) <= eps
                    break; 
                end
            end

            % update Y in the SDP variables
            sdpvars.Y = Y;
        end
    end
    
    % obtain Lyapunov matrix 
        % M, N, P22
        [M,N] = factorize(eye(n)-X*Y',nc1);
        P22 = -M'*Y/N';
        sdpvars.X = X;
        sdpvars.Y = Y; 
        sdpvars.P = [X  M;
                     M' P22];
        
    % useful for the reconstruction
    sdpvars.M = M;
    sdpvars.N = N; 
   
end

%% CONTINUOUS TIME - RECONSTRUCTION

function K = ct_controller_reconstruction(sdpvars,opts)

    plant = sdpvars.plant; 
    nc = size(sdpvars.M,2);
     
    switch lower(opts.reconstruction.method)
        case 'transform'        
            % reconstruct T
            T = [sdpvars.Y*plant.A'+plant.A*sdpvars.Y'  plant.A                                 plant.Bw                sdpvars.Y*plant.Cz' ;
                 plant.A'                               plant.A'*sdpvars.X+sdpvars.X'*plant.A   sdpvars.X'*plant.Bw     plant.Cz'           ;
                 plant.Bw'                              plant.Bw'*sdpvars.X                     -sdpvars.Gammaw         sdpvars.D'          ;
                 plant.Cz*sdpvars.Y'                    plant.Cz                                sdpvars.D               -sdpvars.Gammaz     ];
            
            % reconstruct R
            R = [zeros(plant.nx())      eye(plant.nx())                 zeros(plant.nx(),plant.nw()+plant.nz()); 
                 plant.Bu'              zeros(plant.nu(),plant.nx())    zeros(plant.nu(),plant.nw()) plant.Dzu'];

            % reconstruct S
            S = [eye(plant.nx())                zeros(plant.nx())      zeros(plant.nx(),plant.nw()+plant.nz()); 
                 zeros(plant.ny(),plant.nx())   plant.Cy               plant.Dyw  zeros(plant.ny(),plant.nz())];

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
            
            % do the inverse controller parameter transformation
            Theta = Theta - [sdpvars.X'*plant.A*sdpvars.Y' zeros(plant.nx(),plant.ny()); 
                             zeros(plant.nu(),plant.nx()+plant.ny())                   ]; 
            L = [sdpvars.M            sdpvars.X'*plant.Bu; 
                 zeros(plant.nu(),nc) eye(plant.nu())    ];
            R = [sdpvars.N'           zeros(nc,plant.ny());
                 plant.Cy*sdpvars.Y'  eye(plant.ny())     ];
            Theta = L\Theta/R; 
            
        case 'actual'
            % obtain the required matrices
            P = sdpvars.P;

            % reconstruct R
            Bue = [zeros(plant.nx(),nc)    plant.Bu;
                   eye(nc)                 zeros(nc,plant.nu())]; 
            Dzue = [zeros(plant.nz(),nc)   plant.Dzu];
            R = [Bue'*P	zeros(plant.nu()+nc,plant.nw()) Dzue'];

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
            T = [A0'*P+P'*A0    P'*B0              C0'            ;
                 B0'*P          -sdpvars.Gammaw    D0'            ;
                 C0             D0                 -sdpvars.Gammaz];

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
                    error('LMI solver for feasibility problem (constroller reconstruction) not found. Supported interfaces are ''lmilab'', ''yalmip'', ''cvx'' and ''basiclmi''.'); 
            end
            
        otherwise
            error('Reconstruction method should be ''actual'' (default) or ''transform''.');
            
    end
    
    % cast state-space controller matrices
    Ak = Theta(1:nc,1:nc);
    Bk = Theta(1:nc,(nc+1):end);
    Ck = Theta((nc+1):end,1:nc); 
    Dk = Theta((nc+1):end,(nc+1):end); 
         
    % account for the algebraic loop of the elimination of the impulsive modes
    Mcyu2 = eye(plant.nu())+Dk*plant.Dyu2;
    Ck = Mcyu2\Ck;
    Dk = Mcyu2\Dk;
    Ak = Ak-Bk*plant.Dyu2*Ck;
    Bk = Bk*(eye(plant.ny())-plant.Dyu2*Dk); 
    
    % add the direct feedthrough that eliminates the impulsive modes
    Dk = Dk + plant.Duy; 
    
    % account for the plant's algebraic loop
    Mcyu1 = eye(plant.nu())+Dk*plant.Dyu1;
    Ck = Mcyu1\Ck;
    Dk = Mcyu1\Dk;
    Ak = Ak-Bk*plant.Dyu1*Ck;
    Bk = Bk*(eye(plant.ny())-plant.Dyu1*Dk); 
    
    % transform back to controller for the original plant
    K = hinfcd.dss(Ak,Bk,Ck,Dk,eye(nc),plant.Ts);

end

function Theta = ct_reconstruction_YALMIP(R,S,T,opts)
    disp('Solving with yalmip...');
    % solve for Theta
    t = sdpvar(1,1);
    Theta = sdpvar(size(R,1),size(S,1),'full'); 
    optimize([0.5*(T+T')+R'*Theta*S+S'*Theta'*R <= t*eye(size(T))], t, opts.yalmip);
    Theta = double(Theta); 
end

function Theta = ct_reconstruction_CVX(R,S,T,opts)
    disp('Solving with cvx...');
    % solve for Theta
    cvx_begin SDP
        variable Theta(size(R,1), size(S,1))
        variable t(1,1)
        subject to
            0.5*(T+T')+R'*Theta*S+S'*Theta'*R <= t*eye(size(T))
        eval(['cvx_solver ' opts.cvx.solver]);
        eval(['cvx_precision ' opts.cvx.precision]);
        eval(['cvx_solver_settings(' strjoin(opts.cvx.solver_settings,',') ')']);
    cvx_end 
end

function Theta = ct_reconstruction_LMILAB(R,S,T,opts)
    disp('Solving with lmilab...');
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
    disp('Solving with basiclmi...');
    Theta = basiclmi(T,R,S,'Xmin,Shift');
    %Theta = basiclmi(T,R,S);
end
