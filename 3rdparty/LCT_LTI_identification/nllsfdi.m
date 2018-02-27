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

function [Bn,An,Bls,Als,Bls2,Als2] = nllsfdi(FRF,freq,FRF_W,nh,nl,M_mh,M_ml,iterno,relvar,GN,cORd,varargin)
%
% [Bn,An,Bls,Als,Bls2,Als2]=nllsfdi_v44(FRF,freq,FRF_W,n,M_mh,M_ml,iterno,relvar,GN,cORd,fs)
%
% Non Linear Least Squares FDI (MIMO)
% FRF       : matrix of measured FRF values
%             order the frf's as follows:  [H11 H12 H13 ... H21 H22 H23 .... H31 H32 H33 ...]
%             do this for your own convenience, it is NOT necessary for the program itself because
%             the program does not consider the number of inputs and outputs, it only considers
%             the number of frf's 
% freq      : frequency vector
% FRF_W     : matrix of frequency weighting functions
% nh        : highest order of the denominator polynomial
% nl        : lowest order of the denominator polynomial
% M_mh,M_ml : vector of high and low order of the numerator polynomials
%             the i-th element of M_mh and M_ml corresponds to the i-th FRF
% iterno    : number of iterations
% relvar    : minimum relative deviation of the cost function
% Bn,An     : solution after "iter" iterations
% Bls,Als   : LS-solution
% GN        : if GN==1 : Gauss Newton optimization, otherwise: Levenberg-Marquardt
% cORd      : if 'c', identification of a continuous time model
%             if 'd', identification of a discrete time model          
% fs        : sampling frequency (optional parameter)

% Check if FRF_W is consistent
if isempty(FRF_W) %supply fixed weight of ones
    FRF_W = ones(size(FRF));
else
    assert(all(size(FRF)==size(FRF_W)),'The size of the weights should match the size of the FRF');
end
assert((length(M_mh) == size(FRF,2)) && (length(M_ml) == size(FRF,2)),'size of order of the numerator is inconsistent with the provided FRF');

% make vectors of M_mh and M_ml if they entered as matrices
M_mh=M_mh';
M_ml=M_ml';
M_mh = M_mh(:);
M_ml = M_ml(:);

nino=length(M_ml);
%N : number of FRF measurement points
freq = freq(:);
N = length(freq);

%nrofB : total of numerator coefficients
nrofB=sum(M_mh-M_ml)+nino; 
nrofpar = nrofB+nh-nl;

j=sqrt(-1);

% calculation of the frequency axis
if (cORd == 'c')
   waxis = j*2*pi*freq;
   wscale = median(abs(waxis));
   waxis = waxis/wscale;
elseif (cORd == 'd')
   try 
       fs = varargin{1};
   catch
       error('Sampling frequency was not provided.')
   end
   waxis = exp(j*2*pi*freq/fs);
else
   warning('Time domain is undefined. It is set to continuous time.')
   cORd = 'c';
   waxis = j*2*pi*freq;   
   wscale = median(abs(waxis));
   waxis = waxis/wscale;
end;

if (max(M_mh) > nh)
   warning('The order of the numerator is larger than the order of the denominator.');
end;
if (min(M_mh-M_ml) < 0) 
   error(' \n Some elements of M_lh are larger than the corresponding elements of M_mh. \n');
end;
   

% calculation of the frequency matrices used to avoid the use of polyval
ex=(nh:-1:nl)';
EX=kron(ones(N,1),ex');
Wp = kron(ones(1,nh-nl+1),waxis);
P_fr= Wp.^EX;

maxB = max(M_mh);
minB = min(M_ml);
ex=(maxB:-1:minB)';
EX=kron(ones(N,1),ex');
Wq = kron(ones(1,maxB-minB+1),waxis);
Q_fr = Wq.^EX;

% scaling of the measured FRF with the weighting function
FRF_WD = FRF.*FRF_W;

%
%LLS-schatting = initial estimate for iterative process
% 

fprintf('\n Calculation of the LLS solution. \n')

P_fr = kron(ones(nino,1),P_fr);
P = P_fr.*kron(ones(1,nh-nl+1),FRF_WD(:));

Q = zeros(nino*N,nrofB);
index_count=1;
for i=1:nino
    nr = M_mh(i)-M_ml(i);
    Q(N*(i-1)+1:N*i,index_count:index_count+nr) = kron(ones(1,nr+1),FRF_W(:,i)).*Q_fr(:,maxB-M_mh(i)+1:maxB-M_ml(i)+1);
    index_count = index_count+nr+1;
end; 

A = [real(P(:,2:nh-nl+1)) -real(Q);imag(P(:,2:nh-nl+1)) -imag(Q)];
b = -1*[real(P(:,1)) ; imag(P(:,1))];
y=pinv(A)*b;

[Bls,Als] = BA_construct(y,nh,nl,M_mh,M_ml);

% repeat the least squares estimate with better weighting functions
fprintf('\n Re-calculation of the LLS solution. \n')

invAjw = 1./abs(P_fr(1:N,:)*[1; y(1:nh-nl)]);
FRF_WW = FRF_W.*kron(ones(1,nino),invAjw);

%
%LLS-schatting = initial estimate for iterative process
% 

P = P.*kron(ones(nino,nh-nl+1),invAjw);

index_count=1;
for (i=1:nino)
  Q(N*(i-1)+1:N*i,index_count:index_count+M_mh(i)-M_ml(i)) = kron(ones(1,M_mh(i)-M_ml(i)+1),FRF_WW(:,i)).*Q_fr(:,maxB-M_mh(i)+1:maxB-M_ml(i)+1);
  index_count = index_count + M_mh(i)-M_ml(i)+1;
end; 

A = [real(P(:,2:nh-nl+1)) -real(Q);imag(P(:,2:nh-nl+1)) -imag(Q)];
b = -1*[real(P(:,1)) ; imag(P(:,1))];
y=pinv(A)*b;

% storing the solution in matrices An and Bn
% An is just a row vector
% Bn is a matrix, in which the polynomial coefficients
% are stored

[Bls2,Als2] = BA_construct(y,nh,nl,M_mh,M_ml);

%calculation of cost for Bls2, Als2
invAjw = 1./(P_fr(1:N,:)*[1; y(1:nh-nl)]);

index_count=nh-nl+1;
E=[];
for (i=1:nino)
  nr = M_mh(i)-M_ml(i);
  Ei = (FRF(:,i) - Q_fr(:,maxB-M_mh(i)+1:maxB-M_ml(i)+1)*y(index_count:index_count+nr).*invAjw).*FRF_W(:,i);
  E = [E; Ei];
  index_count = index_count+nr+1;
end;
cost0 = sum(abs(E).^2)/2;

fprintf('\n Start of NLLS iterations \n');
%
% iterative estimation according to the Levenberg-Marquardt algorithm
% 
%
%initialization of parameters
if GN==1
  relax = 0;
else
 relax = 0.01;
end;

iter = 0;
iter0 = 0;
relerror0 = Inf;
relerror = Inf;
y0 = y;
dE_dB = zeros(nino*N,nrofB);

% start of iterative process
% 
while (iter<iterno)&(abs(relerror)>relvar)

  iter = iter + 1;

  Ajw = P_fr(1:N,:)*[1; y(1:nh-nl)]; % Ajw = polyval(An,waxis);

  E = [];
  Bjw = [];
  index_count=nh-nl+1;
  for (i=1:nino)
    nr = M_mh(i)-M_ml(i);  
    Bjw = [Bjw Q_fr(:,maxB-M_mh(i)+1:maxB-M_ml(i)+1)*y(index_count:index_count+nr)]; % Bjw = [Bjw polyval(Bn(i,:),waxis)];   
    E = [E; (FRF(:,i)-Bjw(:,i)./Ajw).*FRF_W(:,i)];
    index_count = index_count+nr+1;
  end;

% scaling of the measured FRF with the weighting function
  FRF_WD = FRF_W./kron(ones(1,nino),Ajw.^2);

  dE_dA=[];
  if (nh-nl~=0)
    dE_dA= P_fr(:,2:end).*kron(ones(1,nh-nl),FRF_WD(:).*Bjw(:));
  end;

  FRF_WD = FRF_W./kron(ones(1,nino),Ajw);

  index_count=1;
  for (i=1:nino)
    nr = M_mh(i)-M_ml(i);    
    dE_dB(N*(i-1)+1:N*i,index_count:index_count+M_mh(i)-M_ml(i))=-kron(ones(1,nr+1),FRF_WD(:,i)).*Q_fr(:,maxB-M_mh(i)+1:maxB-M_ml(i)+1);
    index_count = index_count+nr+1;
  end; 

  A = [real(dE_dA) real(dE_dB);imag(dE_dA) imag(dE_dB)];
  b = [real(E) ; imag(E)];
  JtJ = A'*A;
  Jte = A'*b;
  diagJtJ = diag(JtJ);
   
  A1 = A;
  A2 =  sqrt(relax*diag(diagJtJ+max(diagJtJ)*eps));
  A = [A1 ; A2];
  b = [b; zeros(nrofpar,1)];
  dy = - pinv(A)*b;  
  y0 = y;
  y = y + dy;

% calculation of the cost  
invAjw = 1./(P_fr(1:N,:)*[1; y(1:nh-nl)]);

index_count=nh-nl+1;
E=[];
for (i=1:nino)
  nr = M_mh(i)-M_ml(i);
  Ei = (FRF(:,i) - Q_fr(:,maxB-M_mh(i)+1:maxB-M_ml(i)+1)*y(index_count:index_count+nr).*invAjw).*FRF_W(:,i);
  E = [E; Ei];
  index_count = index_count+nr+1;
end;
cost = sum(abs(E).^2)/2;

relerror = (cost-cost0)/cost0; %relative deviation of the cost
  
if ((cost < cost0)|(GN==1))
   y0 = y; cost0=cost; %updating solution 
   iter0 = iter;relerror0 = relerror;
   relax = relax/2; %lowering the Levenberg-Marquardt factor
else
   y = y0;cost=cost0; %restoring the best result 
   relax = relax*10;
end;    

  fprintf('teller = %g  ITER = %g  cost = %g  rel. change of cost = %g \n',iter,iter0,cost0,relerror0)
 
end;

[Bn,An] = BA_construct(y,nh,nl,M_mh,M_ml);
%rescaling of all models in case of continuous time models
if strcmp(cORd, 'c')
    [no,m]=size(Bn);
    scale = kron(ones(no,1),wscale.^[-m+1:1:0]);
    Bn = Bn.*scale;

    [no,m]=size(An);
    scale = kron(ones(no,1),wscale.^[-m+1:1:0]);
    An = An.*scale;

    [no,m]=size(Bls);
    scale = kron(ones(no,1),wscale.^[-m+1:1:0]);
    Bls = Bls.*scale;

    [no,m]=size(Als);
    scale = kron(ones(no,1),wscale.^[-m+1:1:0]);
    Als = Als.*scale;

    [no,m]=size(Bls2);
    scale = kron(ones(no,1),wscale.^[-m+1:1:0]);
    Bls2 = Bls2.*scale;

    [no,m]=size(Als2);
    scale = kron(ones(no,1),wscale.^[-m+1:1:0]);
    Als2 = Als2.*scale;
end


