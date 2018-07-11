function [Ak,Bk,Ck,Dk] = controller_recon_cL(A,Bu,Cy,Xe,Ye,Ac_hat,Bc_hat,Cc_hat,Dc_hat,param)
% Dimensions
for i = 1:length(param)
    arr(i) = mean(param{i}.range);
end
X = Xe.eval(arr);
Y = Ye.eval(arr);

nx = size(A,1); nu = size(Bu,2); ny = size(Cy,1);                    

% Apply singular value decomposition
[Left,D1,Right] = svd(eye(size(X,1))-X*Y); 
M = Left*sqrt(D1); N = Right*sqrt(D1);

% Transformation
Dk = Dc_hat;
Ck = (Cc_hat-Dk*Cy*X)*pinv(M');
Bk = pinv(N)*(Bc_hat - Y*Bu*Dk);
Ak = pinv(N)*(Ac_hat - N*Bk*Cy*X - Y*Bu*Ck*M' - Y*(A + Bu*Dk*Cy)*X)*pinv(M');
% Extract controller state-space matrices
% Ak = K(1:nx    , 1:nx); Bk = K(1:nx    , nx+1:end);
% Ck = K(nx+1:end, 1:nx); Dk = K(nx+1:end, nx+1:end);
end