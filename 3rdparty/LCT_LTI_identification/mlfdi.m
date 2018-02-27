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

function [Bn,An,Bls,Als,cost0,costls]=mlfdi(X,Y,freq,sX2,sY2,cXY,n,mh,ml,iterno,relvar,cORd,fs)
%
% [Bn,An,Bls,Als,cost0,costls]=mlfdi(X,Y,freq,sX2,sY2,cXY,n,mh,ml,iterno,relvar,cORd,fs)
% 
% Maximum Likelihood FDI 
% X         : input values of the FRF
% Y         : output values of the FRF
% freq      : frequency vector
% sX2       : variance of the input frequency domain noise 
% sY2       : variance of the output frequency domain noise 
% cXY       : covariance of input and output frequency domain noise 
% n         : order of the denominator polynomial
% mh,ml     : high and low order of the numerator polynomial
% iterno    : number of iterations
% relvar    : minimum relative deviation of the cost function
% Bn,An     : solution after "iter" iterations
% Bls,Als   : LS-solution
% cORd      : if 'c', identification of a continuous time model
%             if 'd', identification of a discrete time model          
% fs        : sampling frequency (optional parameter)

j=sqrt(-1);
freq=freq(:);
N=length(freq);
nrofpar = mh-ml+n+1;

% calculation of the frequency axis
if (cORd == 'c')
   waxis = j*2*pi*freq;
elseif (cORd == 'd')
   waxis = exp(j*2*pi*freq/fs);
else
   fprintf(' \n time domain is undefined; it is set to continuous time \n')
   cORd = 'c';
   waxis = j*2*pi*freq;
end;

%
% LS-schatting = initial estimate for iterative process
% 

fprintf(' \n calculation of the LS solution \n')

% matrices P and Q are used to form the complex set of equations
P = kron(ones(1,n+1),waxis).^kron(ones(N,1),(n:-1:0));
Q = kron(ones(1,mh-ml+1),waxis).^kron(ones(N,1),(mh:-1:ml));

%calculation of the matrices A and b : A * y = b
PP = kron(ones(1,n+1),Y).*P;
QQ = kron(ones(1,mh-ml+1),X).*Q;

A = [real(PP(:,2:n+1)) -real(QQ);imag(PP(:,2:n+1)) -imag(QQ)];
b = -1*[real(PP(:,1)) ; imag(PP(:,1))];
y=pinv(A)*b;

% store the solution in matrices An and Bn
An = [1 y(1:n)'];
Bn = [zeros(1,n-mh) y(n+1:nrofpar)' zeros(1,ml)];


% least squares solution is called Bls, and Als
Als = An;
Bls = Bn;

%calculation of cost for Bls, Als
cost0 = mlfdi_res(Bls,Als,freq,X,Y,sX2,sY2,cXY,cORd,fs);
fprintf('ML cost function for least squares solution  = %g   \n',cost0);
costls = cost0;

fprintf('\n start of iteration \n');
%
% iterative estimation according to the Levenberg-Marquardt algorithm
% 
%
%initialization of parameters
relax = 1;
skip = 0;
iter = 0;
relerror0 = inf;
relerror = inf;
y0 = y;

% start of iterative process
% 
while (iter<iterno)&(relerror>relvar)

  iter = iter + 1;

%
%calculation of new set of equations
%

  Num = Q*y(n+1:n+1+mh-ml);
  Den = P*[1;y(1:n)];
  SE = sqrt( sX2.*(abs(Num).^2) + sY2.*(abs(Den).^2) - 2*real(cXY.*Den.*conj(Num)) );
  E = Num.*X - Den.*Y;

  A = [];

  for (i=2:n+1)
     WW = -Y.*P(:,i)./SE - E./(SE.^3).*( sY2.*real(Den.*conj(P(:,i)))- real(cXY.*P(:,i).*conj(Num)));
     A = [A WW];
  end;   


  for (i=1:mh-ml+1)
     WW = X.*Q(:,i)./SE - E./(SE.^3).*( sX2.*real(Num.*conj(Q(:,i))) - real(cXY.*Den.*conj(Q(:,i))) );
     A = [A WW];
  end;   

  J = [real(A);imag(A)];
  f = [real(E)./SE ; imag(E)./SE];
   
  JtJ = J'*J;
  Jte = J'*f;
  diagJtJ = diag(JtJ);
   
  A1 = J;
  A2 =  sqrt(relax*diag(diagJtJ+max(diagJtJ)*eps));
  A = [A1 ; A2];
  b = [f; zeros(nrofpar,1)];
  dy = -pinv(A)*b;   
  y0 = y;
  y = y + dy;

% storing the solution in matrices An and Bn
An = [ 1 y(1:n)'];
Bn = [zeros(1,n-mh) y(n+1:nrofpar)' zeros(1,ml)];

% calculation of the cost 
  cost = mlfdi_res(Bn,An,freq,X,Y,sX2,sY2,cXY,cORd,fs);
 
  relerror = abs(cost-cost0)/cost0; %relative deviation of the cost
  
  if (cost < cost0)
     y0 = y; cost0=cost; %updating solution 
     iter0 = iter;relerror0 = relerror;
     relax = relax/2; %lowering the Levenberg-Marquardt factor
  else
     y = y0;cost=cost0; %restoring the best result 
     An = [ 1 y(1:n)'];
     Bn = [zeros(1,n-mh) y(n+1:nrofpar)' zeros(1,ml)];
     relax = relax*10;
  end;    

  fprintf('teller = %g  ITER = %g  cost = %g  relative eror = %g \n',iter,iter0,cost0,relerror0)
 
end;