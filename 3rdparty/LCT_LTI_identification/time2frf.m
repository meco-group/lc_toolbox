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

function [X,Y,freq,sX2,sY2,cXY,FRF,sCR] = time2frf(x,y,fs,nl,nh,nrofsapp)
%
% Version 2: addition of FRF and sCR (20-02-2002)
%
% [X,Y,freq,sX2,sY2,cXY,FRF,sCR] = time2frf(x,y,fs,nl,nh,nrofsapp)
%
% x, y     : input and output data vector measured using periodic broad-band
%            excitation 
% fs       : sampling frequency
% nl, nh   : number of lowest and highest frequency point of the excitated 
%            frequency band. If DC is excited: nl = 0.
% nrofsapp : number of sample points per period
% X, Y     : input and output frf-values : FRF = Y./X
% sX2, sY2 : variance of real and imaginary parts of noise on X and Y values
% cXY      : covariance between real and imaginary parts of noise on X and Y
% FRF      : frequency response function
% sCR      : cramer-rao variance on FRF
%

x=x(:);
y=y(:);


freq= (nl:1:nh)'/(nrofsapp/fs);

nrofp=floor(length(x)/nrofsapp);
nr=length(freq);

INP=[];
OUT=[];
for (i=1:nrofp)
	INP=[INP fft(x(1+(i-1)*nrofsapp:i*nrofsapp))];
	OUT=[OUT fft(y(1+(i-1)*nrofsapp:i*nrofsapp))];
end;

INP=INP(nl+1:1:nh+1,:);
OUT=OUT(nl+1:1:nh+1,:);




X = INP(:,1);
Y = OUT(:,1);
FRF=Y./X;
sX2 = zeros(nr,1);
sY2 = zeros(nr,1);
cXY = zeros(nr,1);
sCR = zeros(nr,1);



if (nrofp>1)
    X=(mean(INP'))';
    Y=(mean(OUT'))';
    sX2=((std(INP'))'.^2)/2/nrofp;
    sY2=((std(OUT'))'.^2)/2/nrofp;
    
    for (i=1:nr)
        Q=cov(INP(i,:),OUT(i,:));
        cXY(i,1) = Q(1,2)/2/nrofp;
    end;
    FRF=Y./X;
    sCR=2*abs(FRF).*(sX2./(abs(X)).^2+sY2./(abs(Y)).^2-2*real(cXY./(conj(X).*Y)));
end;

