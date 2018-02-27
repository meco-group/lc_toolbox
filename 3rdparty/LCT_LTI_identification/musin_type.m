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

function [x,X,freq,Xt,freqt] = musin_type(fs,fmin,fmax,nrofsamp,initial,cp,tp,Bn,An)
%
% [x,X,freq,Xt,freqt] = musin_type(fs,fmin,fmax,nrofsamp,initial,cp,tp,Bn,An)
%
% fs : sample frequency in Hz
% fmin : minimal frequency
% fmax : maximal frequency
% nrofsamp : number of samples per period
% x : multi sine with flat amplitude spectrum (time domain signal)
% X : spectrum of multi sine in frequency band [fmin,fmax]
% freq : frequency axis in Hz
% Xt : total spectrum of multi sine  
% freqt : total frequency axis in Hz
% initial: initial guess of phases (optional): 
%		's' is schroeder phases (defaul
%		'r' is random phases between -pi and pi
% cp : 'c' = compressed multi sine (default), 'n' = not compressed multi sine
% tp : type: 'o' = odd, 'O' = odd-odd, 'O2' = special odd-odd, 'f' = full multi sine (also default)
% Bn,An	: continuous time transfer function model describing the desired amplitude spectrum: | Bn / An |
%			: default: flat spectrum ...
%

if nargin<5, initial = 's'; end
if initial~='r', initial = 's'; end
if nargin<6, cp = 'c'; end
if nargin<7, tp = 'full'; end
if nargin<8, Bn = 1; An = 1; end

df = fs/nrofsamp;
n_max=round(fmax/df)+1;
n_min=ceil(fmin/df)+1;
w_s=j*2*pi*df*(0:1:nrofsamp/2-1)';
w_s=w_s(1:n_max);


Mag = abs(polyval(Bn,w_s)./polyval(An,w_s));

if (tp ~= 'o' & tp~='O') %construction of full multi sine flat spectrum
   X0=zeros(n_max,1);
   X0(n_min:1:n_max) = Mag(n_min:1:n_max);
   ff = 1;
end;

if tp == 'o'   %odd multi sine flat spectrum
   vc = (1:2:n_max-1)';
   i = find(abs(vc-(n_min-1)) == min(abs(vc-(n_min-1))));
   i = max(i);
   ind = vc(i)+1;
   X0=zeros(n_max,1); 
   X0(ind:2:n_max) = Mag(ind:2:n_max);%ones(length(ind:2:n_max),1);
   ff = 2;   
end;

if tp == 'O'
   %odd-odd multi sine flat spectrum
   vc = (1:4:n_max-1)';
   i = find(abs(vc-(n_min-1)) == min(abs(vc-(n_min-1))));
   i = max(i);
   ind = vc(i)+1;
   X0=zeros(n_max,1); %construction of the flat spectrum
   X0(ind:4:n_max) = Mag(ind:4:n_max);%ones(length(ind:4:n_max),1);
   ff = 4;   
end;



if strcmp(tp,'O2')
   %odd-odd multi sine flat spectrum
   vc1 = (1:8:n_max-1)';
   vc2 = (3:8:n_max-1)';
   vc = sort([vc1;vc2]);
   i = find(abs(vc-(n_min-1)) == min(abs(vc-(n_min-1))));
   i = max(i);
   X0=zeros(n_max,1); %construction of the flat spectrum
   X0(ind:4:n_max) = Mag(ind:4:n_max);%ones(length(ind:4:n_max),1);
   ind = vc(i)+1;
   if max(vc(i)==vc1)      
      X0(ind:8:n_max) = Mag(ind:8:n_max);
      X0(ind+2:8:n_max) = Mag(ind+2:8:n_max);
   end;
   if max(vc(i)==vc2)
      X0(ind:8:n_max) = Mag(ind:8:n_max);
      X0(ind+6:8:n_max) = Mag(ind+6:8:n_max);
   end;
   ff = 4;   
end;

p=2; % starting value of p
X=X0;
if (cp ~= 'n');
   while p<600      % Minimizes the l-2p norm
      X=msinl2pi(p,X,nrofsamp,[],0,[],10,1e-4,initial);
      % Increment the value of p
      p=ceil(p*2);
   end      % stops when p>=600
else
   if (initial=='r') 
      X = randph(X); % If X=real then randph
   else
      X=schroed(X); % If X=real then schroefrt
   end;
end;

%figure
x = four2ti(X,nrofsamp); % calculation of the time domain multi-sine
X=time2fo(x,nrofsamp); %recalculation of the spectrum
freq=fs*(0:1:nrofsamp/2-1)'/nrofsamp;
%subplot(2,1,1);
%plot(freq(2:ff:nrofsamp/2),20*log10(abs(X(2:ff:nrofsamp/2))));xlabel('frequency in Hz'); ylabel('amplitude in dB');
%title('total frequency characteristic of the designed multi-sine');
%subplot(2,1,2);
time=(0:1/fs:(nrofsamp-1)/fs)';
%plot(time,x); xlabel('time in s'); title('multi-sine');
freqs=freq(n_min:n_max);
Xs = X(n_min:n_max);
Xt = X;
freqt = freq;
X = Xs;
freq = freqs;
%subplot(111)