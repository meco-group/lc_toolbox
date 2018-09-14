function [ y ] = simLPVc( varargin )

required = {'model','data'};
optional = {'x0'};

calls = varargin(1:2:end);
assert(mod(nargin,2) == 0,'Odd number of arguments.');
assert(iscellstr(calls),'The names of your name/value pairs should be strings.');
assert(~any(~ismember(required,calls)),['These arguments are required: ' strjoin(required)]);
for i = 1:nargin/2
    if strcmp(varargin{2*i-1},'model'); assert(~exist('model','var'),'You cannot define model twice.'); model = varargin{2*i}; end;
    if strcmp(varargin{2*i-1},'data'); assert(~exist('data','var'),'You cannot define data twice.'); data = varargin{2*i}; end;
    if strcmp(varargin{2*i-1},'x0'); assert(~exist('x0','var'),'You cannot define x0 twice.'); x0 = varargin{2*i}; end;
    if ~ismember(varargin{2*i-1},[required, optional]); warning(['I could not identify the option ''' varargin{2*i-1} ''' and I am going to ignore it.']); end;
end

if ~exist('x0','var')
    x0 = zeros(size(model_new.A.coeff_tensor,2),1);
end

basisEval = model.A.basis.list_eval(globalData.sch); 
cfg = coder.config('mex');
cfg.DynamicMemoryAllocation = 'AllVariableSizeArrays';
codegen -config cfg evalnstackGlob -args { zeros(size(model.A.coeff_tensor)), zeros(size(model.B.coeff_tensor)), zeros(size(model.C.coeff_tensor)), zeros(size(model.D.coeff_tensor)), zeros(size(permute(basisEval,[2 3 1]))) }                                                                                           
[Ap,Bp,Cp,Dp] = stackBases_mex(model.A.coeff_tensor,model.B.coeff_tensor,model.C.coeff_tensor,model.D.coeff_tensor,permute(basisEval,[2 3 1]));
codegen -config cfg simLPV -args { zeros(size(Ap)), zeros(size(Bp)), zeros(size(Cp)), zeros(size(Dp)), zeros(size(globalData.u)), zeros(size(model.A.coeff_tensor,2),1) }                                                                        

y = simLPV_mex(Ap,Bp,Cp,Dp,data.u,x0);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ yc ] = simLPV( A,B,C,D,u,x0 )

% n  -> model order, nu -> nr of inputs, ny -> nr of outputs, N  -> nr of data samples 
% size(A) = N n n, size(B) = N n nu, size(C) = N ny n, size(D) = N ny nu, size(u) = N nu 1, size(x0) =  n 1                                                                     
x = x0;
yc = zeros(size(u,1),size(C,2),1);
for t = 1:size(u,1)
    yc(t,:,:) = reshape(C(t,:,:),[size(C,2) size(C,3)])*x + reshape(D(t,:,:),[size(D,2) size(D,3)])*reshape(u(t,:,:),[size(u,2) size(u,3)]);
    x = reshape(A(t,:,:),[size(A,2) size(A,3)])*x + reshape(B(t,:,:),[size(B,2) size(B,3)])*reshape(u(t,:,:),[size(u,2) size(u,3)]); 
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ As, Bs, Cs, Ds ] = stackBases( A, B, C, D, basisEval)
% suppose size(basisEval) = [nb 1 N]
As = zeros(size(basisEval,3),size(A,2),size(A,3));
Bs = zeros(size(basisEval,3),size(B,2),size(B,3));
Cs = zeros(size(basisEval,3),size(C,2),size(C,3));
Ds = zeros(size(basisEval,3),size(D,2),size(D,3));

for t = 1:size(basisEval,3)
    for k = 1:size(basisEval,1)
        As(t,:,:) = As(t,:,:) + A(k,:,:)*basisEval(k,:,t);
        Bs(t,:,:) = Bs(t,:,:) + B(k,:,:)*basisEval(k,:,t);
        Cs(t,:,:) = Cs(t,:,:) + C(k,:,:)*basisEval(k,:,t);
        Ds(t,:,:) = Ds(t,:,:) + D(k,:,:)*basisEval(k,:,t);
    end 
end

end
