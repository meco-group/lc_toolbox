function [ As, Bs, Cs, Ds ] = evalnstackGlob( A, B, C, D, basisEval)
% suppose size(basisEval) = [nb 1 N]
As = zeros(size(basisEval,3),size(A,2),size(A,3));
Bs = zeros(size(basisEval,3),size(B,2),size(B,3));
Cs = zeros(size(basisEval,3),size(C,2),size(C,3));
Ds = zeros(size(basisEval,3),size(D,2),size(D,3));

for t = 1:size(basisEval,3)
    for i = 1:size(basisEval,1)
        As(t,:,:) = As(t,:,:) + A(i,:,:)*basisEval(i,:,t);
        Bs(t,:,:) = Bs(t,:,:) + B(i,:,:)*basisEval(i,:,t);
        Cs(t,:,:) = Cs(t,:,:) + C(i,:,:)*basisEval(i,:,t);
        Ds(t,:,:) = Ds(t,:,:) + D(i,:,:)*basisEval(i,:,t);
    end
    
end

end

