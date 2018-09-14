function [ y ] = simLPV( A,B,C,D,u,x0 )

% n  -> model order, nu -> nr of inputs, ny -> nr of outputs, N  -> nr of data samples 
% size(A) = N n n, size(B) = N n nu, size(C) = N ny n, size(D) = N ny nu, size(u) = N nu 1, size(x0) =  n 1
% for codegen: 
% cfg = coder.config('mex');
% cfg.DynamicMemoryAllocation = 'AllVariableSizeArrays';
% codegen -config cfg simLPV -args { zeros(size(A)), zeros(size(B)), zeros(size(C)), zeros(size(D)), zeros(size(u)), zeros(size(x0)) }                                                                        
x = x0;
y = zeros(size(u,1),size(C,2),1);
for t = 1:size(u,1)
    y(t,:,:) = reshape(C(t,:,:),[size(C,2) size(C,3)])*x + reshape(D(t,:,:),[size(D,2) size(D,3)])*reshape(u(t,:,:),[size(u,2) size(u,3)]);
    x = reshape(A(t,:,:),[size(A,2) size(A,3)])*x + reshape(B(t,:,:),[size(B,2) size(B,3)])*reshape(u(t,:,:),[size(u,2) size(u,3)]); 
end

end