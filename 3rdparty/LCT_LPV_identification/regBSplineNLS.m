function varargout = regBSplineNLS( varargin )

required = {'model','schParam','gamma'};
optional = {'griddedFRFs','globalData','FRFweights','alpha','paramActivity','scaleCones','reweighting','relWtol','minRelKdiff','maxLambda','epsilon','maxIter','minDSnorm'};

calls = varargin(1:2:end);
assert(mod(nargin,2) == 0,'Odd number of arguments.');
assert(iscellstr(calls),'The names of your name/value pairs should be strings.');
assert(~any(~ismember(required,calls)),['These arguments are required: ' strjoin(required)]);
for i = 1:nargin/2
    if strcmp(varargin{2*i-1},'model'); assert(~exist('model','var'),'You cannot define model twice.'); model = varargin{2*i}; end;
    if strcmp(varargin{2*i-1},'schParam'); assert(~exist('schParam','var'),'You cannot define schParam twice.'); schParam = varargin{2*i}; end;
    if strcmp(varargin{2*i-1},'gamma'); assert(~exist('gamma','var'),'You cannot define gamma twice.'); gamma = varargin{2*i}; end;   
    if strcmp(varargin{2*i-1},'griddedFRFs'); assert(~exist('griddedFRFs','var'),'You cannot define griddedFRFs twice.'); griddedFRFs = varargin{2*i}; end;
    if strcmp(varargin{2*i-1},'globalData'); assert(~exist('globalData','var'),'You cannot define globalData twice.'); globalData = varargin{2*i}; end;
    if strcmp(varargin{2*i-1},'FRFweights'); assert(~exist('FRFweights','var'),'You cannot define FRFweights twice.'); FRFweights = varargin{2*i}; end;
    if strcmp(varargin{2*i-1},'alpha'); assert(~exist('alpha','var'),'You cannot define alpha twice.'); alpha = varargin{2*i}; end;
    if strcmp(varargin{2*i-1},'paramActivity'); assert(~exist('paramActivity','var'),'You cannot define paramActivity twice.'); paramActivity = varargin{2*i}; end;
    if strcmp(varargin{2*i-1},'scaleCones'); assert(~exist('scaleCones','var'),'You cannot define scaleCones twice.'); scaleCones = varargin{2*i}; end;
    if strcmp(varargin{2*i-1},'reweighting'); assert(~exist('reweighting','var'),'You cannot define reweighting twice.'); reweighting = varargin{2*i}; end;
    if strcmp(varargin{2*i-1},'relWtol'); assert(~exist('relWtol','var'),'You cannot define relWtol twice.'); relWtol = varargin{2*i}; end;
    if strcmp(varargin{2*i-1},'minRelKdiff'); assert(~exist('minRelKdiff','var'),'You cannot define minRelKdiff twice.'); minRelKdiff = varargin{2*i}; end;
    if strcmp(varargin{2*i-1},'minDSnorm'); assert(~exist('minDSnorm','var'),'You cannot define minDSnorm twice.'); minDSnorm = varargin{2*i}; end;
    if strcmp(varargin{2*i-1},'maxLambda'); assert(~exist('maxLambda','var'),'You cannot define maxLambda twice.'); maxLambda = varargin{2*i}; end;
    if strcmp(varargin{2*i-1},'epsilon'); assert(~exist('epsilon','var'),'You cannot define epsilon twice.'); epsilon = varargin{2*i}; end;
    if strcmp(varargin{2*i-1},'maxIter'); assert(~exist('maxIter','var'),'You cannot define epsilon twice.'); maxIter = varargin{2*i}; end;
    if ~ismember(varargin{2*i-1},[required, optional]); warning(['I could not identify the option ''' varargin{2*i-1} ''' and I am going to ignore it.']); end;
end
assert(isa(model,'LPVDSSmod'),'Unsupported data type.')

meco_binaries('cpp_splines','develop')
import splines.*

if exist('griddedFRFs','var') 
    p_local = model.A.basis.list_eval(griddedFRFs.params_{2})';
    cfg = coder.config('mex');
    cfg.DynamicMemoryAllocation = 'AllVariableSizeArrays';
    codegen -config cfg evalnstackLoc -args { zeros(size(model.A.coeff_tensor)), zeros(size(model.B.coeff_tensor)), zeros(size(model.C.coeff_tensor)), zeros(size(model.D.coeff_tensor)), zeros(size(permute(p_local,[1 3 2]))) }                                                                                           
    freqLines = 1;
    OmegaConc = [];
    FRFmConc  = [];  
    for k = 1:numel(griddedFRFs.params_{2})
        OmegaConc = [OmegaConc; 2*pi*griddedFRFs.grid_{k}.Frequency];
        FRFmConc = [FRFmConc; vec(griddedFRFs.grid_{k}.frdata)];
        freqLines = [freqLines; (freqLines(end) + numel(griddedFRFs.grid_{k}.Frequency))];
    end
    [As,Bs,Cs,Ds] = evalnstackLoc_mex(model.A.coeff_tensor,model.B.coeff_tensor,model.C.coeff_tensor,model.D.coeff_tensor,permute(p_local,[1 3 2]));
    codegen -config cfg getFRFjacobian_LCT_rearr -args { zeros(size(p_local,2),size(model.A.coeff_tensor,2),size(model.A.coeff_tensor,3)), zeros(size(p_local,2),size(model.B.coeff_tensor,2),size(model.B.coeff_tensor,3)), zeros(size(p_local,2),size(model.C.coeff_tensor,2),size(model.C.coeff_tensor,3)), zeros(size(p_local,2),size(model.D.coeff_tensor,2),size(model.D.coeff_tensor,3)),0, zeros(numel(griddedFRFs.grid_)+1,1), zeros(numel(griddedFRFs.grid_)*numel(griddedFRFs.grid_{1}.frdata),1), zeros(size(p_local)) }                                                                                           
    if exist('FRFweights','var')
        Wfreq = sqrt(1./sumsqr(abs(FRFmConc).*FRFweights))*diag(repmat(FRFweights,2,1));
    else
       % Wfreq = diag((1/sqrt(sumsqr(abs(FRFmConc))))*ones(1,2*numel(FRFmConc)));
        FRFweights = ones(size(FRFmConc))./(abs(FRFmConc)); % comment out for the crane example
        Wfreq = sqrt(1./sumsqr(abs(FRFmConc).*FRFweights))*diag(repmat(FRFweights,2,1)); % comment out for the crane example
    end
    if exist('alpha','var')
        Wfreq = sqrt(alpha)*Wfreq; 
    else
        if exist('globalData','var')
            Wfreq = sqrt(0.5)*Wfreq;
        end
    end
    if any(any(Wfreq > 0))
        errFRF = Wfreq*getFRFerror_LCT(freqLines,OmegaConc,FRFmConc,griddedFRFs.params_{2},schParam,model);
        Jf = Wfreq*getFRFjacobian_LCT_rearr_mex(As,Bs,Cs,Ds,model.Ts,freqLines,OmegaConc,p_local);
    else
        errFRF = [];
        Jf = [];
    end
else
    errFRF = [];
    Jf = [];
end

if exist('globalData','var') % to be adapted (for the time being assumed that globalData is a structure with fields (tensors) u ([N nu 1]), sch ([N 1 1]) and y([N ny 1]))     

    basisEval = model.A.basis.list_eval(globalData.sch);
    cfg = coder.config('mex');
    cfg.DynamicMemoryAllocation = 'AllVariableSizeArrays';
    codegen -config cfg evalnstackGlob -args { zeros(size(model.A.coeff_tensor)), zeros(size(model.B.coeff_tensor)), zeros(size(model.C.coeff_tensor)), zeros(size(model.D.coeff_tensor)), zeros(size(permute(basisEval,[2 3 1]))) }                                                                                           
    [Ap,Bp,Cp,Dp] = evalnstackGlob_mex(model.A.coeff_tensor,model.B.coeff_tensor,model.C.coeff_tensor,model.D.coeff_tensor,permute(basisEval,[2 3 1]));
    codegen -config cfg simLPV -args { zeros(size(Ap)), zeros(size(Bp)), zeros(size(Cp)), zeros(size(Dp)), zeros(size(globalData.u)), zeros(size(model.A.coeff_tensor,2),1) }                                                                        
    codegen -config cfg getJacobianFast_LCT -args { zeros(size(Ap)), zeros(size(Bp)), zeros(size(Cp)), zeros(size(globalData.u)), zeros(size(permute(basisEval,[2 3 1]))), zeros(size(model.A.coeff_tensor,2),1) }                                                                        
    Wt = 1/sqrt(sumsqr(globalData.y));
    if exist('alpha','var')
        Wt = sqrt(1-alpha)*Wt; 
    else
        if exist('griddedFRFs','var')
            Wt = sqrt(0.5)*Wt;
        end
    end
    if Wt > 0
        errTD = Wt*(simLPV_mex(Ap,Bp,Cp,Dp,globalData.u,zeros(size(model.A.coeff_tensor,2),1)) - globalData.y);
        Jt =  Wt*getJacobianFast_LCT_mex(Ap,Bp,Cp,globalData.u,permute(basisEval,[2 3 1]),zeros(size(model.A.coeff_tensor,2),1));
    else 
        errTD = [];
        Jt = [];
    end
else
    errTD = [];
    Jt = []; 
end

err = [errFRF; errTD];
J = [Jf; Jt];
K = err'*err;

opti = OptiSplineYalmip();
ops = sdpsettings('solver','mosek','savesolveroutput',1,'verbose',0);
ops.mosek.MSK_DPAR_INTPNT_CO_TOL_DFEAS = 1e-12;
ops.mosek.MSK_DPAR_INTPNT_CO_TOL_PFEAS = 1e-12;
opti.solver('yalmip',struct('yalmip_options',ops,'use_optimize',true));
     
T = [model.A model.B; model.C model.D];
S_sizes = [dimension(basis(T)) size(T)];
Delta_S = opti.variable(prod(S_sizes), 1);
S_k = opti.parameter(prod(S_sizes), 1);
sk = opti.parameter(T.basis.dimension-T.basis.degree-1,1);

jacScalTensor = ones(size(T.coeff_tensor));
if exist('paramActivity','var')
    for i = 1:size(T.coeff_tensor,1)
       jacScalTensor(i,:,:) = paramActivity;    
    end
end  
activeIndices = vec(jacScalTensor);
J = J(:,find(activeIndices));

T_Delta_S = Function(basis(T),MTensor(Delta_S,S_sizes));
T_Sk = Function(basis(T),MTensor(S_k,S_sizes));
T_symbolic = T_Sk + T_Delta_S;
sExact = diff(coeff_tensor(T.derivative(T.basis.degree)));

if exist('scaleCones','var')
    if strcmp(scaleCones,'yes')
        coneScaling = squeeze(rms(coeff_tensor(T.derivative(T.basis.degree))));
    elseif strcmp(scaleCones,'no')
        coneScaling = ones(size(sExact,2),size(sExact,3));
    else
        disp('Unsupported entry for scaleCones.')
    end  
else
    coneScaling = ones(size(sExact,2),size(sExact,3));
end

indices = find(coneScaling > 1e-12);
coneScaling_robust = ones(size(sExact,2),size(sExact,3));
coneScaling_robust(indices) = coneScaling(indices);
coneScaling = 1./(coneScaling_robust);

if exist('paramActivity','var')
    coneScaling = paramActivity.*coneScaling;
end

coneScalTensor = ones(size(sExact,1),size(sExact,2),size(sExact,3));
for i = 1:size(coneScalTensor,1)
   coneScalTensor(i,:,:) = coneScaling;   
end
sExact = coneScalTensor.*diff(coeff_tensor(T.derivative(T.basis.degree)));
ds = opti.variable(T.basis.dimension-T.basis.degree-1,1);
s = ones(size(sExact,1),1);
for i = 1:numel(s)
        s(i) = norm(vec(sExact(i,:,:)));
    end
for i = 1:numel(s)
   s(i) 
end 
part = 100*abs(s)./sum(abs(s));

if ~exist('relWtol','var')
    relWtol = 1e-4;
end
if ~exist('minRelKdiff','var')
    minRelKdiff = 1e-8;
end
if ~exist('minDSnorm','var')
    minDSnorm = 1e-3;
end
if ~exist('maxLambda','var')
    maxLambda = 1e8;
end
if ~exist('maxIter','var')
    maxIter = 1e3;
end

relWdiff = Inf;
iteration = 0;
disp([num2str(iteration) ' - cost: ' num2str(K) ' - s: ' num2str(s') ' - stable: ' num2str(isstable(model))]);
disp('Starting the optimization...')

while any(relWdiff > relWtol) && (iteration < maxIter)
    
    lambda = -1; 
    redExit = false;
    relKdiff = Inf;
    cones = coneScalTensor.*diff(coeff_tensor(T_symbolic.derivative(T.basis.degree)));
    w = ones(1,size(sExact,1));
    dSnorm1 = Inf;
        
    while (not(redExit) && (dSnorm1 > minDSnorm))

        if lambda == -1
            Kreg = 0;
            for i = 1:numel(s)
                Kreg = Kreg + gamma*w(i)*norm(vec(sExact(i,:,:)));
            end
            [Uj,Sj,Vj] = svd(J,'econ');
            lambda = sqrt(Sj(1,1)); 
            clear Uj Sj Vj    
        end

        if (isnan(K) || isinf(lambda)), break; end  
        greenExit = false;

        opti.set_value(T_Sk, T);
        opti.set_value(sk, s);
        for i = 1:numel(ds)
            Conei = matrix(cones(i,:,:));
            Conei = Conei(:);
            opti.subject_to(norm(Conei) <= (sk(i) + ds(i)));
        end

        while (not(greenExit) && not(redExit))              
            
            M = 2*(J'*J + (lambda^2)*eye(size(J,2))); 
            opti.minimize((2*J'*err)'*Delta_S + (1/2)* Delta_S'*M*Delta_S + gamma*w*ds); 
            sol = opti.solve();
            Tnew = sol.value(T_symbolic);
            model_new = SSmod(Tnew(1:model.A.size(1),1:model.A.size(2)),Tnew(1:model.B.size(1),model.A.size(2)+1:model.A.size(2)+model.B.size(2)),Tnew(model.A.size(1)+1:model.A.size(1)+model.C.size(1),1:model.C.size(2)),Tnew(model.A.size(1)+1:model.A.size(1)+model.D.size(1),model.C.size(2)+1:model.C.size(2)+model.D.size(2)),schParam,model.Ts);

            if exist('griddedFRFs','var') && (any(any(Wfreq > 0)))
                errFRF_new = Wfreq*getFRFerror_LCT(freqLines,OmegaConc,FRFmConc,griddedFRFs.params_{2},schParam,model_new);
            else
                errFRF_new = [];
            end
            if exist('globalData','var') && (Wt > 0)
                [Ap,Bp,Cp,Dp] = evalnstackGlob_mex(model_new.A.coeff_tensor,model_new.B.coeff_tensor,model_new.C.coeff_tensor,model_new.D.coeff_tensor,permute(basisEval,[2 3 1]));
                errTD_new = Wt*(simLPV_mex(Ap,Bp,Cp,Dp,globalData.u,zeros(size(model_new.A.coeff_tensor,2),1)) - globalData.y);
            else
                errTD_new = [];
            end
            err_new = [errFRF_new; errTD_new];
            Knew = err_new'*err_new;
            Kreg_new = 0;
            sExact_new  = sol.value(cones);
            for i = 1:numel(ds)
                Kreg_new = Kreg_new + gamma*w(i)*norm(vec(sExact_new(i,:,:)));
            end
            if ((Knew + Kreg_new) >= (K + Kreg))
                lambda = lambda*sqrt(10); 
                disp('Searching... ')
                if (lambda >= maxLambda)
                    redExit = true;
                end
            else
                lambda = lambda/2;
                greenExit = true;
            end  

        end

        if greenExit 
            relKdiff = abs((Knew+Kreg_new) - (K+Kreg))./abs(K+Kreg);
            K = Knew;
            Kreg = Kreg_new;
            err = err_new;
            T = Tnew;
            model = model_new; 
            sExact = sExact_new;
            
            if exist('griddedFRFs','var') && (any(any(Wfreq > 0)))
                [As,Bs,Cs,Ds] = evalnstackLoc_mex(model.A.coeff_tensor,model.B.coeff_tensor,model.C.coeff_tensor,model.D.coeff_tensor,permute(p_local,[1 3 2]));
                Jf = Wfreq*getFRFjacobian_LCT_rearr_mex(As,Bs,Cs,Ds,model.Ts,freqLines,OmegaConc,p_local);
            else
                Jf = [];
            end
            if exist('globalData','var') && (Wt > 0)
                [Ap,Bp,Cp,Dp] = evalnstackGlob_mex(model.A.coeff_tensor,model.B.coeff_tensor,model.C.coeff_tensor,model.D.coeff_tensor,permute(basisEval,[2 3 1]));
                Jt =  Wt*getJacobianFast_LCT_mex(Ap,Bp,Cp,globalData.u,permute(basisEval,[2 3 1]),zeros(size(model.A.coeff_tensor,2),1));
            else
                Jt = []; 
            end
            J = [Jf; Jt];
            J = J(:,find(activeIndices));

            s = s + sol.value(ds);
            iteration = iteration + 1;
            dSnorm1 = norm(sol.value(Delta_S));
            disp([num2str(iteration) ' - K: ' num2str(K)  ' - Kreg: ' num2str(Kreg) ' s: ' num2str(s') ' - relKdiff: ' num2str(relKdiff), ' - dSnorm: ' num2str(dSnorm1) ' - stable: ' num2str(isstable(model))]);     
            
        elseif redExit
                disp(['Lambda too big: ' num2str(lambda)]); 
        end
            
    end
    if ~exist('reweighting','var') || strcmp(reweighting,'no')
        break
    end
    if exist('scaleCones','var')
        if strcmp(scaleCones,'yes')
            coneScaling = squeeze(rms(coeff_tensor(T.derivative(T.basis.degree))));
        elseif strcmp(scaleCones,'no')
            coneScaling = ones(size(sExact,2),size(sExact,3));
        else
            disp('Unsupported entry for scaleCones.')
        end  
    else
        coneScaling = ones(size(sExact,2),size(sExact,3));
    end
    indices = find(coneScaling > 1e-12);
    coneScaling_robust = ones(size(sExact,2),size(sExact,3));
    coneScaling_robust(indices) = coneScaling(indices);
    coneScaling = 1./coneScaling_robust;
    if exist('paramActivity','var')
        coneScaling = paramActivity.*coneScaling;
    end
    coneScalTensor = ones(size(sExact,1),size(sExact,2),size(sExact,3));
    for i = 1:size(coneScalTensor,1)
       coneScalTensor(i,:,:) = coneScaling; 
    end
    sExact = coneScalTensor.*diff(coeff_tensor(T.derivative(T.basis.degree)));
    
    w_new = ones(1,numel(s));
    for i = 1:numel(s)
        w_new(i) = 1/(norm(vec(sExact(i,:,:))) + epsilon);
    end 
    relWdiff = abs(w_new - w)./abs(w);
    w = w_new;
    
    for i = 1:numel(s)
        s(i) = norm(vec(sExact(i,:,:)));
    end
        
    if any(relWdiff > relWtol)
        disp('-------------------------------------------------------------------------------------------------------------------------------------');
        fprintf(1,'\n');
        disp(['Reweighting - ||relDw||: ' num2str(norm(relWdiff,1))]);  
        w
    else
        fprintf(1,'\n');
        disp(['No reweighting - ||relDw||: ' num2str(norm(relWdiff,1))]);  
        break
    end
    
end             
disp('End of the optimization.')

varargout{1} = model;
varargout{2} = s;
