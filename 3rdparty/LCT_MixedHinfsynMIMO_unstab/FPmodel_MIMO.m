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

% This file creates and tunes a first principles model
% 090811 aliasing

clear all
close all
clc

% load G090813
% Gfrf = G;
% clear G

%% plant
num_modes = 6;

% parameters
rho = 7.7e3;
b = 0.024;%0.02
h = 0.001;%0.002
L = 0.5;
m = rho*b*h;
E = 2e11;%2e11
I = b*h^3/12;
Ts = 0.001;

% number of elements
n = 100;

% length of elements
l = L/n;

% diagonals K_e
kd_e = [  6*l, -12,   6*l,  12,    0,    0,     0
          0,   2*l^2, -6*l, 4*l^2, 6*l,  0,     0
          0,   0,     -6*l, 12,    -6*l, -12,   0
          0,   0,     0,    4*l^2, -6*l, 2*l^2, 6*l ];

% diagonals M_e
md_e = [  -13*l, 54,     22*l,  156,   0,     0,      0
          0,     -3*l^2, 13*l,  4*l^2, 22*l,  0,      0
          0,     0,      -22*l, 156,   13*l,  54,     0
          0,     0,      0,     4*l^2, -22*l, -3*l^2, -13*l ];

% form K and M
kd = zeros(2*n+2, 7);
md = zeros(2*n+2, 7);

for e = 1:n
    kd(2*e-1:2*e+2,1:7) = kd(2*e-1:2*e+2,1:7) + kd_e;
    md(2*e-1:2*e+2,1:7) = md(2*e-1:2*e+2,1:7) + md_e;
end

K = spdiags(kd, -3:3, 2*n+2, 2*n+2)*E*I/l^3;
M = spdiags(md, -3:3, 2*n+2, 2*n+2)*rho*b*h*l/420;

% add springs connected with environment
spring = 1.00000035;
K(1, 1) = spring*K(1, 1);
K(2*n+1, 2*n+1) = spring*K(2*n+1, 2*n+1);

% solve eigenvalue problem
OPTS.disp = 0;
[ Phi, Lambda ] = eig(full(K), full(M)); % [ Phi, Lambda ] = eigs(K, M, n, 'SM', OPTS);
[ Lambda, S ] = sort(diag(Lambda));
Phi = Phi(:,S);

% tune eigenfrequencies
omegard = sqrt(Lambda)/1.24; % correction for difference model and plant
omegard(5) = 185*2*pi;
omegard(7) = 438*2*pi;

% input and output matrices
actoffset = 13;
sensoffset = 14;
zoffset = 0;

B = zeros(2*(n+1), 2);
B([ 1+2*actoffset ], 1) = ones(1,1);
B([  2*(n-actoffset)+1 ], 2) = ones(1,1);
Cz = zeros([ 2, 2*(n+1) ]);
Cz(1, n+1+2*zoffset) = 1;
Cz(2, n+1+2*zoffset+1) = 1;
Cy = zeros([ 2, 2*(n+1) ]);
Cy(1, [ 1+2*sensoffset, 2*(n-sensoffset)+1 ]) = ones(1,2)/2;
Cy(2, [ 1+2*sensoffset, 2*(n-sensoffset)+1 ]) = [1 -1]/L;
% coordinate transformation to modal form
Mm = Phi.'*M*Phi;
Km = Phi.'*K*Phi;
Bm = Phi.'*B;
sensz = 500;
sensy = 500;
Cmz = sensz*Cz*Phi;
Cmy = sensy*Cy*Phi;

% proportional damping
damping = 0.01;
zeta = ones(size(omegard))*damping;
zeta(1:2) = ones(2,1)*0.1;
zeta(3:4) = 0.01;

% state space model
Ai = zeros(2*num_modes);
Bi = zeros(2*num_modes, 2);
Ciz = zeros(2, 2*num_modes);
Ciy = zeros(2, 2*num_modes);

for m = 1:num_modes
    Ai([2*m-1, 2*m],[2*m-1, 2*m]) = [ 0, 1; -omegard(m)^2, -2*zeta(m)*omegard(m) ];
    Bi([2*m-1, 2*m],:) = [ zeros(1,2); Bm(m,:) ];
    Ciz(:,[2*m-1, 2*m]) = [ Cmz(:,m), zeros(2,1) ];
    Ciy(:,[2*m-1, 2*m]) = [ Cmy(:,m), zeros(2,1) ];
    
    % tune contribution modes
    if m == 1 || m == 2 || m == 3
        Ciy(:,[2*m-1, 2*m]) = 0.9*Ciy(:,[2*m-1, 2*m]);
    end
    
    if m == 3
        Bi([2*m-1, 2*m],:) = 0.95*Bi([2*m-1, 2*m],:);
    end
    
    if m == 5 || m == 7
        Ciz(:,[2*m-1, 2*m]) = 0.36*Ciz(:,[2*m-1, 2*m]);
        Ciy(:,[2*m-1, 2*m]) = 0.64*Ciy(:,[2*m-1, 2*m]);
    end
end

Gzy = ss(Ai, Bi, [ Ciz; Ciy ], 0);
% Gzy.InputName = {'u1', 'u2'};
% Gzy.OutputName = {'z','zt', 'y1','y2'};

beam_model_MIMO3=Gzy;
save 'beam_model_MIMO3' beam_model_MIMO3

% % collocated control
% Gzy = minreal(balreal(Gzy(:,1)));
% 
% % construct delay
% s = tf('s');
% th = 0%0.05*1e-3; % should be 0 <= th < Ts
% subdelay = exp(-s*th);
% subdelay = ss(subdelay);
% onedelay = ss(0, 1, 1, 0, Ts);
% Gzyd = c2d(Gzy*subdelay, Ts, 'zoh')*onedelay; % add some onedelay's
% % clear Gzy
% Gfrf.Ts = Ts;
% 
% disp('Check if delay is included in the a,b,c,d state space description!');
% 
% %% optimal control
% om_bw = 50*2*pi;
% W2 = c2d(tf([ 1/(om_bw/5), 1 ], [ 1/(om_bw*5), 1 ])*...
%     tf([ 1/(om_bw/6), 1 ], [ 1, 0 ])*...
%     tf(1, [ 1/(om_bw*5), 1 ]), Ts, 'tustin');
% W2 = W2/abs(freqresp(W2*Gzyd(2), om_bw));
% % W1 = tf(1,1,Ts);
% % % Kopt = c2d(-ncfsyn(d2c(Gzyd(2), 'tustin'), d2c(W2, 'tustin')), Ts, 'tustin');
% % Kopt = ncfsyn(minreal(Gzyd(2)), W2);
% Gs = balreal(minreal(W2*Gzyd(2)));
% systemnames = 'Gs';
% inputvar = '[ ry; ru; u ]';
% outputvar = '[ Gs; ru+u; ry-Gs ]';
% input_to_Gs = '[ ru+u ]';
% Ps = sysic;
% Ks_opt = hinfsyn(Ps, 1, 1, 'tolgam', 1e-5, 'display', 'on');
% Kopt = Ks_opt*W2;
% 
% T = feedback(minreal(Gzyd*Kopt), 1, 1, 1);
% step(T); return
% 
% %% plot
% figure;
% Gfrf.InputNames = {'u'};
% Gzyd.InputNames = {'u'};
% 
% bode(Gfrf, 'b-', Gzyd, 'g--', (1:500)*2*pi);
% Ch = get(gcf, 'Children');
% % pause
% legend(Ch(end), 'Non-parametric model', 'Finite element model');
% set(gcf, 'PaperPositionMode', 'auto');
% % print -depsc2 FEMmodel.eps