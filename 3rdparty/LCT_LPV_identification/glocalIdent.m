function varargout = glocalIdent( varargin )

required = {'model','schParam'};
optional = {'griddedFRFs','globalData','FRFweights','alpha','paramActivity','maxLambda','minGradNorm','maxIter'};

calls = varargin(1:2:end);
assert(mod(nargin,2) == 0,'Odd number of arguments.');
assert(iscellstr(calls),'The names of your name/value pairs should be strings.');
assert(~any(~ismember(required,calls)),['These arguments are required: ' strjoin(required)]);
for i = 1:nargin/2
    if strcmp(varargin{2*i-1},'model'); assert(~exist('model','var'),'You cannot define model twice.'); model = varargin{2*i}; end;
    if strcmp(varargin{2*i-1},'schParam'); assert(~exist('schParam','var'),'You cannot define schParam twice.'); schParam = varargin{2*i}; end;
    if strcmp(varargin{2*i-1},'griddedFRFs'); assert(~exist('griddedFRFs','var'),'You cannot define griddedFRFs twice.'); griddedFRFs = varargin{2*i}; end;
    if strcmp(varargin{2*i-1},'globalData'); assert(~exist('globalData','var'),'You cannot define globalData twice.'); globalData = varargin{2*i}; end;
    if strcmp(varargin{2*i-1},'FRFweights'); assert(~exist('FRFweights','var'),'You cannot define FRFweights twice.'); FRFweights = varargin{2*i}; end;
    if strcmp(varargin{2*i-1},'alpha'); assert(~exist('alpha','var'),'You cannot define alpha twice.'); alpha = varargin{2*i}; end;
    if strcmp(varargin{2*i-1},'paramActivity'); assert(~exist('paramActivity','var'),'You cannot define paramActivity twice.'); paramActivity = varargin{2*i}; end;
    if strcmp(varargin{2*i-1},'maxLambda'); assert(~exist('maxLambda','var'),'You cannot define maxLambda twice.'); maxLambda = varargin{2*i}; end;
    if strcmp(varargin{2*i-1},'minGradNorm'); assert(~exist('minGradNorm','var'),'You cannot define minGradNorm twice.'); minGradNorm = varargin{2*i}; end;
    if strcmp(varargin{2*i-1},'maxIter'); assert(~exist('maxIter','var'),'You cannot define maxIter twice.'); maxIter = varargin{2*i}; end;
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
    codegen -config cfg getFRFjacobian_LCT -args { zeros(size(p_local,2),size(model.A.coeff_tensor,2),size(model.A.coeff_tensor,3)), zeros(size(p_local,2),size(model.B.coeff_tensor,2),size(model.B.coeff_tensor,3)), zeros(size(p_local,2),size(model.C.coeff_tensor,2),size(model.C.coeff_tensor,3)), zeros(size(p_local,2),size(model.D.coeff_tensor,2),size(model.D.coeff_tensor,3)),0, zeros(numel(griddedFRFs.grid_)+1,1), zeros(numel(griddedFRFs.grid_)*numel(griddedFRFs.grid_{1}.frdata),1), zeros(size(p_local)) }                                                                                           
    if exist('FRFweights','var')
        Wfreq = sqrt(1./sumsqr(abs(FRFmConc).*FRFweights))*diag(repmat(FRFweights,2,1)); 
    else
       % Wfreq = diag((1/sqrt(sumsqr(abs(FRFmConc))))*ones(1,2*numel(FRFmConc))); %
       FRFweights = ones(size(FRFmConc))./(abs(FRFmConc));
       Wfreq = sqrt(1./sumsqr(abs(FRFmConc).*FRFweights))*diag(repmat(FRFweights,2,1)); 
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
        Jf = Wfreq*getFRFjacobian_LCT_mex(As,Bs,Cs,Ds,model.Ts,freqLines,OmegaConc,p_local);
    else
        errFRF = [];
        Jf = [];
    end   
    
else
    disp('No local data used.')
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
    codegen -config cfg getJacobianFast_LCT_arr -args { zeros(size(Ap)), zeros(size(Bp)), zeros(size(Cp)), zeros(size(globalData.u)), zeros(size(permute(basisEval,[2 3 1]))), zeros(size(model.A.coeff_tensor,2),1) }                                                                        
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
        Jt =  Wt*getJacobianFast_LCT_arr_mex(Ap,Bp,Cp,globalData.u,permute(basisEval,[2 3 1]),zeros(size(model.A.coeff_tensor,2),1));
    else
        errTD = [];
        Jt = [];
    end 
    %Jcasadi = Wt*getJacobian_casadi(model.A.coeff_tensor,model.B.coeff_tensor,model.C.coeff_tensor,model.D.coeff_tensor,globalData.u,basisEval,zeros(size(model.A.coeff_tensor,2),1)); 
else
    errTD = [];
    Jt = [];
    disp('No global data used.')
end
err = [errFRF; errTD];
J = [Jf; Jt];
K = err'*err;

activeParameters = ones(1,numel(model.A.coeff_tensor)+numel(model.B.coeff_tensor)+numel(model.C.coeff_tensor)+numel(model.D.coeff_tensor));
if exist('paramActivity','var')
    Aact = paramActivity(1:size(model.A,1),1:size(model.A,2));
    Bact = paramActivity(1:size(model.B,1),size(model.A,2)+1:size(model.A,2)+size(model.B,2));
    Cact = paramActivity(size(model.A,1)+1:size(model.A,1)+size(model.C,1),1:size(model.C,2));
    Dact = paramActivity(size(model.B,1)+1:size(model.B,1)+size(model.D,1),size(model.C,2)+1:size(model.C,2)+size(model.D,2));
    for i = 1:size(model.A.coeff_tensor,1)
       activeParameters(1,(i-1)*size(model.A,1)*size(model.A,2)+1:i*size(model.A,1)*size(model.A,2)) = Aact(:); 
       activeParameters(1,size(model.A.coeff_tensor,1)*size(model.A,1)*size(model.A,2)+(i-1)*size(model.B,1)*size(model.B,2)+1:size(model.A.coeff_tensor,1)*size(model.A,1)*size(model.A,2)+i*size(model.B,1)*size(model.B,2)) = Bact(:);
       activeParameters(1,size(model.A.coeff_tensor,1)*(size(model.A,1)*size(model.A,2)+size(model.B,1)*size(model.B,2))+(i-1)*size(model.C,1)*size(model.C,2)+1:size(model.A.coeff_tensor,1)*(size(model.A,1)*size(model.A,2)+size(model.B,1)*size(model.B,2))+i*size(model.C,1)*size(model.C,2)) = Cact(:); 
       activeParameters(1,size(model.D.coeff_tensor,1)*(size(model.A,1)*size(model.A,2)+size(model.B,1)*size(model.B,2)+size(model.C,1)*size(model.C,2))+(i-1)*size(model.D,1)*size(model.D,2)+1:size(model.D.coeff_tensor,1)*(size(model.A,1)*size(model.A,2)+size(model.B,1)*size(model.B,2)+size(model.C,1)*size(model.C,2))+i*size(model.D,1)*size(model.D,2)) = Dact(:); 
    end
end      
J = J(:,find(activeParameters));

if ~exist('maxLambda','var')
    maxLambda = 1e8;
end
if ~exist('maxIter','var')
    maxIter = 1e3;
end

lambda = -1;
iteration = 0;
gradNorm = norm(2*J'*err);
disp([num2str(iteration) ' - cost: ' num2str(K)  ' - gradNorm: ' num2str(gradNorm) ' - stable: ' num2str(isstable(model))]);
disp('Starting the optimization...')

while ((gradNorm > minGradNorm) && (lambda < maxLambda) && (iteration < maxIter)) 

    if lambda == -1
        [Uj,Sj,Vj] = svd(J,'econ');
        lambda = sqrt(Sj(1,1)); 
        clear Uj Sj Vj    
    end

    if (isnan(K) || isinf(lambda)), break; end  
    [J0,scalingJ0] = normalise(J);
    nowexit = false;

    while not(nowexit) 
        
        q = - (J0'*J0 + (lambda^2)*eye(size(J0,2))) \ (J0'*err);
        dtheta = q./scalingJ0';
        
        dcoeff_tensor_A = permute(reshape(dtheta(1:numel(model.A.coeff_tensor)),size(model.A.coeff_tensor,2),size(model.A.coeff_tensor,3),size(model.A.coeff_tensor,1)),[3 1 2]);
        dcoeff_tensor_B = permute(reshape(dtheta((numel(model.A.coeff_tensor)+1):(numel(model.A.coeff_tensor)+numel(model.B.coeff_tensor))),size(model.B.coeff_tensor,2),size(model.B.coeff_tensor,3),size(model.B.coeff_tensor,1)),[3 1 2]);
        dcoeff_tensor_C = permute(reshape(dtheta((numel(model.A.coeff_tensor)+numel(model.B.coeff_tensor)+1):(numel(model.A.coeff_tensor)+numel(model.B.coeff_tensor)+numel(model.C.coeff_tensor))),size(model.C.coeff_tensor,2),size(model.C.coeff_tensor,3),size(model.C.coeff_tensor,1)),[3 1 2]);
        dcoeff_tensor_D = permute(reshape(dtheta((numel(model.A.coeff_tensor)+numel(model.B.coeff_tensor)+numel(model.C.coeff_tensor)+1):(numel(model.A.coeff_tensor)+numel(model.B.coeff_tensor)+numel(model.C.coeff_tensor)+numel(model.D.coeff_tensor))),size(model.D.coeff_tensor,2),size(model.D.coeff_tensor,3),size(model.D.coeff_tensor,1)),[3 1 2]);
        A_new = splines.Function(model.A.tensor_basis,model.A.coeff_tensor + dcoeff_tensor_A);
        B_new = splines.Function(model.B.tensor_basis,model.B.coeff_tensor + dcoeff_tensor_B);
        C_new = splines.Function(model.C.tensor_basis,model.C.coeff_tensor + dcoeff_tensor_C);
        D_new = splines.Function(model.B.tensor_basis,model.D.coeff_tensor + dcoeff_tensor_D);
        model_new = SSmod(A_new,B_new,C_new,D_new,schParam,model.Ts);

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
            
        if (Knew >= K)
            lambda = lambda*sqrt(10); 
        else
            lambda = lambda/2;
            nowexit = true;
        end  

    end

    iteration = iteration + 1;
    K = Knew;
    err = err_new;
    model = model_new; 
            
    if exist('griddedFRFs','var') && (any(any(Wfreq > 0))) 
        [As,Bs,Cs,Ds] = evalnstackLoc_mex(model.A.coeff_tensor,model.B.coeff_tensor,model.C.coeff_tensor,model.D.coeff_tensor,permute(p_local,[1 3 2]));
        Jf = Wfreq*getFRFjacobian_LCT_mex(As,Bs,Cs,Ds,model.Ts,freqLines,OmegaConc,p_local);
    else
        Jf = [];
    end
    if exist('globalData','var') && (Wt > 0)
        [Ap,Bp,Cp,Dp] = evalnstackGlob_mex(model.A.coeff_tensor,model.B.coeff_tensor,model.C.coeff_tensor,model.D.coeff_tensor,permute(basisEval,[2 3 1]));
        Jt =  Wt*getJacobianFast_LCT_arr_mex(Ap,Bp,Cp,globalData.u,permute(basisEval,[2 3 1]),zeros(size(model.A.coeff_tensor,2),1));
      %  Jt =  Wt*getJacobianFast_LCT_mex(Ap,Bp,Cp,globalData.u,permute(basisEval,[2 3 1]),zeros(size(model.A.coeff_tensor,2),1));
    else
        Jt = []; 
    end
    J = [Jf; Jt];
    J = J(:,find(activeParameters));
    
    gradNorm = norm(2*J'*err);
    disp([num2str(iteration) ' - K: ' num2str(K) ' - grad: ' num2str(gradNorm) ' - lambda: ' num2str(lambda) ' - stable: ' num2str(isstable(model))]);     
     
end             
disp('End of the optimization.')
varargout{1} = model;
varargout{2} = err;

