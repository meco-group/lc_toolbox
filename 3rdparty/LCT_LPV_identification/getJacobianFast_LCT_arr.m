function [Jacobian] = getJacobianFast_LCT_arr(As, Bs, Cs, u, basisEval, x0)

% size(As) = N n n, size(Bs) = N n nu 
% size(Cs) = N ny n, size(Ds) = N ny nu
% size(u) = N nu 1
% size(basisEval) = nb 1 N 
na  = size(As,2);
nu = size(Bs,3);
ny = size(Cs,2);
nb = size(basisEval,1);
N  = size(u,1);
x = zeros(N+1,na,1);
x(1,:,1) = x0;
dx = zeros(N+1,na,nb*(na*na+na*nu+ny*na+ny*nu));
Jacobian = zeros(ny*N,nb*(na*na+na*nu+ny*na+ny*nu));

for t = 1:N      
    for i = 1:nb
        % A
        for k = 1:na
            for j = 1:na
                Jacobian(((t-1)*ny+1):(t*ny),(i-1)*na*na+(k-1)*na+j) = reshape(Cs(t,:,:),[ny na])*reshape(dx(t,:,(i-1)*na*na+(k-1)*na+j),[na 1]);
                dx(t+1,:,(i-1)*na*na+(k-1)*na+j) = reshape(basisEval(i,:,t),[1 1])*Ijk([na na],[j k])*reshape(x(t,:,:),[na 1]) + reshape(As(t,:,:),[na na])*reshape(dx(t,:,(i-1)*na*na+(k-1)*na+j),[na 1]);  
            end       
        end
        % B
        for k = 1:nu
            for j = 1:na
                Jacobian(((t-1)*ny+1):(t*ny),nb*na*na+(i-1)*na*nu+(k-1)*na+j) = reshape(Cs(t,:,:),[ny na])*reshape(dx(t,:,nb*na*na+(i-1)*na*nu+(k-1)*na+j),[na 1]);
                dx(t+1,:,nb*na*na+(i-1)*na*nu+(k-1)*na+j) = reshape(basisEval(i,:,t),[1 1])*Ijk([na nu],[j k])*reshape(u(t,:,:),[size(u,2) size(u,3)]) + reshape(As(t,:,:),[na na])*reshape(dx(t,:,nb*na*na+(i-1)*na*nu+(k-1)*na+j),[na 1]); 
            end       
        end
        % C
        for k = 1:na
            for j = 1:ny
                Jacobian(((t-1)*ny+1):(t*ny),nb*(na*na+na*nu)+(i-1)*ny*na+(k-1)*ny+j) = reshape(basisEval(i,:,t),[1 1])*Ijk([ny na],[j k])*reshape(x(t,:,:),[na 1]);
            end       
        end
        % D
        for k = 1:nu
            for j = 1:ny
                Jacobian(((t-1)*ny+1):(t*ny),nb*(na*na+na*nu+ny*na)+(i-1)*ny*nu+(k-1)*ny+j) = reshape(basisEval(i,:,t),[1 1])*Ijk([ny nu],[j k])*reshape(u(t,:,:),[size(u,2) size(u,3)]);
            end       
        end
    end 
    x(t+1,:,:) = reshape(As(t,:,:),[na na])*reshape(x(t,:,:),[na 1]) + reshape(Bs(t,:,:),[na nu])*reshape(u(t,:,:),[size(u,2) size(u,3)]);
end

end

