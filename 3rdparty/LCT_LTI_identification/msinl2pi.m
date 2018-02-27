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

function X=msinl2pi(p,X,Nmax,Fa,W,H,iterno,relvar,initial)

%MSINL2P This function minimizes the L2p-norm of a multisine.

% Use X=msinl2p(p,X,Nmax,Fa,W,H,iterno,relerr,initial) with

%  Input arguments :

%  - p  : Value of p in L2p

%  - X  : Matrix with input spectra (1=DC)

%  - Nmax : Max. number of points in time domain

%  - Fa : Set of the additional harmonic numbers (DC = 1)

%         (auto power spectrum not restricted)

%  - W  : Weighting row vector, e.g., W=[1 2];

%         shortcuts :  W=0 -> same crest factors (default)

%                      W=1 -> same amplitudes

%  - H  : Transfer function matrix 

%         (no.row=no.freq.;no.col.=1 (SISO version))

%  - iterno : Max. number of iterations

%  - relvar : Max. relative variation of L2p-norm

%  Output arguments :

%  - X  : Matrices with optimal input spectra

% initial: initial guess of phases (optional): 

%		's' is schroeder phases (default)

%		'r' is random phases between -pi and pi

%

% Example :      X=ones(32,1);X(1)=0; % X(1) -> DC

%                p=2;

%                while p<200

%                   X=msinl2p(p,X,2048,[],0,[],10,1e-6);

%                   p=ceil(p*2);

%                end

%

% See also SCHROED, CREST1, CREST2.

 %

%       P.A.N. Guillaume - version 1 / 24 November 1990

%       Copyright (c) 1990 by dept. ELEC, V.U.B.

%

if nargin<8, relvar=1e-6; end      % Default rel. dev. cost f.

if nargin<7, iterno=10; end        % 10 iterations by default

if nargin<6, H=[]; end             % Only input cf-minimization

if nargin<5, W=0; end, W=W(:)';    % W is a row vector (or scalar)

if nargin<4, Fa=[]; end, Fa=Fa(:); % Fa is a column vector

if nargin<3, Nmax=2048; end

if nargin<9, initial = 's'; end

if initial~='r', initial = 's'; end



%       Fe : Set of the effective harmonic numbers (index 1 = DC)

[rowno,colno]=size(X);  % Fe = all indices of non-zero X's entries - Fa

if colno==1, dummy=X; else dummy=max(abs(X)')'; end

Fe=cumsum(ones(size(dummy)));

Fe=Fe(dummy~=0);

clear dummy

for I=1:length(Fa), Fe(Fe==Fa(I))=[]; end  % Elimination of Fa indices

Feno=length(Fe); % Number of effective frequencies

Fano=length(Fa); % Number of additional frequencies

Ft=[Fe;Fa];      % Set of all the considered spectral lines

Ftno=Feno+Fano;  % Total number of considered spectral lines

kmax=max([max(Fe),max(Fa)])-1;  % Largest harmonic number of X

if length(X(:,1))<=kmax, X(kmax+1,1)=0; end

Ndummy=2*p*kmax+1; % Number of points needed for identity between 

     % L2p (summation) and l2p (integration) criterions

N=2; while Ndummy>=N, N=N*2; end % Select N=2^x > Ndummy

N=min([N,Nmax]);                 % If N>Nmax then N:=Nmax



if all(imag(X)==0), 

  if (initial=='r') 

    X = randph(X); % If X=real then randph

  else

    X=schroed(X); % If X=real then schroed

  end;

end;

    

if ~isempty(H)

   H=H(1:kmax+1);% Correct range setting of H

   Y=zeros(size(X));

   Y=H.*X;           % Output spectra

   if W==0, W=[effval(X,Fe),effval(Y,Fe)].^(-1); end

   % Same crest factors

   if W==1, W=[effval(X,Fe),effval(X,Fe)].^(-1); end

   % Same extreme value

   X=X*W(1);Y=Y*W(2);

   x=four2ti(X,N);               % Input signals

   y=four2ti(Y,N);               % Output signals

else

   W=1/effval(X,Fe);

   X=X*W;        % Scaling so that max|x(t)|~1

   x=four2ti(X,N);               % Input signals

end

if isempty(H)

   cost=lnorm(x(:),2*p);        % Calculation of the L2p-norm

else

   cost=lnorm([x(:);y(:)],2*p);

end

    % (L2p=costfunction)

X0=X;cost0=cost;iter0=0;        % Initialization of the best results

    % (Best = lowest cost function)

indexmat=(Ft-1)*ones(size(Ft))';      % Initialization of the (k-l)-matrix

indexmin=indexmat-indexmat'+N;  % and the (k+l)-matrix (Offset=N)

indexplus=indexmat+indexmat'+N;

phkphlmin=indexmin(:);          % Transformation of the matrices

phkphlplus=indexplus(:);        % into a vector

if ~isempty(Fa)

   phkalmin=indexmin(:,Feno+1:Ftno);

   phkalplus=indexplus(:,Feno+1:Ftno);

   akalmin=indexmin(Feno+1:Ftno,Feno+1:Ftno);

   akalplus=indexplus(Feno+1:Ftno,Feno+1:Ftno);

   phkalmin=phkalmin(:);        % Transformation of the matrices

   phkalplus=phkalplus(:);      % into a vector

   akalmin=akalmin(:);

   akalplus=akalplus(:);   

end

relax=0.01;skip=0;              % Levenberg-Marquardt parameters

iter=0;                         % Iteration counter

relerror=inf;                   % Relative deviation of the cost function

%/// MAIN ITERATION LOOP ////////////////////////////////////////////////

while (iter<iterno)&(relerror>relvar)  % Stop criterions

   iter=iter+1;

   X2p2=fft(x.^(2*p-2));               % FFT of x^(2p-2)

   X2p2=[flipud(conj(X2p2(2:N)));X2p2];

   dummyt=X(Ft);

   %  Calculation of JphJph and JphE

   Qmin=conj(dummyt*dummyt');

   Qplus=conj(dummyt*dummyt.');

   Qmin(:)=Qmin(:).*X2p2(phkphlmin);

   Qplus(:)=Qplus(:).*X2p2(phkphlplus);

   JphJph=p*real(Qmin-Qplus);

   JphE=sum(imag(Qmin+Qplus)')';

   if ~isempty(Fa)

      dummya=exp(j*angle(X(Fa)));

      %  Calculation of JphJa and JaE

      Qmin=conj(dummyt*dummya');

      Qplus=conj(dummyt*dummya.');

      Qmin(:)=Qmin(:).*X2p2(phkalmin);

      Qplus(:)=Qplus(:).*X2p2(phkalplus);

      JphJa=p*imag(Qmin+Qplus);

      JaE=sum(real(Qmin+Qplus))';

      %  Calculation of JaJa

      Qmin=conj(dummya*dummya');

      Qplus=conj(dummya*dummya.');

   Qmin(:)=Qmin(:).*X2p2(akalmin);

      Qplus(:)=Qplus(:).*X2p2(akalplus);

      JaJa=p*real(Qmin+Qplus);

   end

   if ~isempty(H)

      Y=(H.*X)*W(2)/W(1);

      y=four2ti(Y,N);

      Y2p2=fft(y.^(2*p-2));            % FFT of y^(2p-2)

      Y2p2=[flipud(conj(Y2p2(2:N)));Y2p2];

      dummyt=Y(Ft);

      %  Calculation of JphJph and JphE

      Qmin=conj(dummyt*dummyt');

      Qplus=conj(dummyt*dummyt.');

      Qmin(:)=Qmin(:).*Y2p2(phkphlmin);

      Qplus(:)=Qplus(:).*Y2p2(phkphlplus);

      JphJph=JphJph+p*real(Qmin-Qplus);

      JphE=JphE+sum(imag(Qmin+Qplus)')';

      if ~isempty(Fa)

  dummya=exp(j*angle(Y(Fa)));

  %  Calculation of JphJa and JaE

  Qmin=conj(dummyt*dummya');

  Qplus=conj(dummyt*dummya.');

         Qmin(:)=Qmin(:).*Y2p2(phkalmin);

  Qplus(:)=Qplus(:).*Y2p2(phkalplus);

  JphJa=JphJa+p*imag(Qmin+Qplus);

  JaE=JaE+sum(real(Qmin+Qplus))';

  %  Calculation of JaJa

  Qmin=conj(dummya*dummya');

  Qplus=conj(dummya*dummya.');

      Qmin(:)=Qmin(:).*Y2p2(akalmin);

  Qplus(:)=Qplus(:).*Y2p2(akalplus);

  JaJa=JaJa+p*real(Qmin+Qplus);

      end

   end

   if isempty(Fa)

      A=JphJph;b=JphE;

   else

      A=[JphJph,JphJa;JphJa',JaJa];b=[JphE;JaE];

   end

   diagA=diag(A);       

   A=A+relax*diag(diagA+max(diagA)*eps); % Adding diagonal matrix to A

      % (Levenberg-Marquardt)

   Delta=-A\b;

   X(Fe)=abs(X(Fe)).*exp(j*(angle(X(Fe))+Delta(1:Feno))); % Updating of X

   X(Fa)=(abs(X(Fa))+Delta(Ftno+1:Ftno+Fano)).*exp(j*(angle(X(Fa))+Delta(Feno+1:Ftno)));

   if ~isempty(H),Y=(H.*X)*W(2)/W(1);y=four2ti(Y,N);end% Updating of Y and y

   x=four2ti(X,N);                                                         

                                                                           

                                                                           

                                                % Updating of x

   if isempty(H)

      cost=lnorm(x(:),2*p);            % Calculation of the L2p-norm

   else

      cost=lnorm([x(:);y(:)],2*p);

   end

   relerror=abs(cost-cost0)/cost0;     % Relative deviation of the cost f.

   if cost<cost0

      X0=X;cost0=cost;  % Updating X0 (the best multi-sine spectra)

      iter0=iter;

      relax=relax/2;    % Lowering the Levenberg-Marquardt factor

   else

      X=X0;cost=cost0;  % Restoring the best results

      x=four2ti(X,N); 

      if ~isempty(H), Y=(H.*X)*W(2)/W(1); y=four2ti(Y,N); end

      relax=relax*10;   % Augmenting the Levenberg-Marquardt factor

   % (Bad convergence)

   end

  % fprintf('P=%g/CF=%g/ITER=%g',p,crest2(X,Nmax,Fe),iter)

  % fprintf('/NORM=%g/ITER0=%g/NORM0=%g',cost,iter0,cost0)

  % fprintf('/LM=%g\n',relax)

end

X=X0/W(1);  % Unscaling the result



