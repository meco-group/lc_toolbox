%% Velocity controller design for mechatronic drive train
%  Laurens Jacobs, 18 Oct 2019

% initialize
clear; close all; clc; 
addpath(genpath('lc_toolbox'));

%% identification

% nonparametric data
load('frequencyRange.mat'); 
Ga = [];
for i=1:6
    load(['FilterSpeedFRF161019run' num2str(i) '.mat']);
    G{i} = FRDmod(FRF,f,'FrequencyUnit','Hz');
end
Gav = 1/4*(G{1}+G{2}+G{3}+G{4}); % 5 and 6 are saturating
figure; bode(G{1:4}); hold on; bode(Gav,'k'); title('Measurements & average'); 

% fit
s = tf('s'); 
Gtofit = Gav*s; 
FRFW = squeeze(ones(size(Gtofit.ResponseData))); 
FRFW(1:100) = 0;        % don't fit to 10 Hz
FRFW(300:800) = 1.5;    % append more weight to the AR/R pair
FRFW(1900:end) = 0;     % don't fit from 190 Hz
settings = struct('denh',4,'denl',0,'numh',2,'numl',0,'FRFW',FRFW);
Gfit = 1/s*param_ident('data',Gtofit,'method','nllsfdi','settings',settings);
figure; bode(Gav); hold on; bode(Gfit); title('Model'); 

% calculate relative error of fit (multiplicative uncertainty estimate)
for i=1:4
    relerr{i} = std((Gfit-G{i}))/std(Gfit); 
end
figure; bodemag(relerr{:});

%% control design

% control configuration
P = IOSystem(Gfit);         % plant
K = IOSystem(1,1);          % controller

r = Signal(1);              % velocity reference
u = P.in;                   % current command
y = P.out;                  % measured speed
e = r-y;                    % tracking error
conn = [K.in == e; 
        K.out == u];    
CL = IOSystem(P,K,conn);    % define the control configuration

% define channels
S = Channel(e/r,'Sensitivity'); 
T = Channel(y/r,'Complementary sensitivity');
U = Channel(u/r,'Input sensitivity');

% define specifications
WS = Weight.LF(10,2,-80); 
MS = Weight.DC(6);
WU = Weight.HF(20,2,-80);

obj =  WS*S;
cstr = [MS*S <= 1; 
        WU*U <= 1];
    
% solve
[~,~,info] = CL.solve(obj,cstr,K);

% plot
% showall(info); % plots the channels that were in the optimization problem
% figure; bode(CL(T)); title('Closed loop'); % plot the closed-loop
% figure; bode(K); title('Controller'); % plot the controller

% postprocess
[z,p,k] = zpkdata(std(K)); 
z = sort(z{:}); p = sort(p{:});
z(1) = []; p(1) = [];
p(1) = 0;
k = k/prod(abs(p(abs(p)>30*2*pi))); p(abs(p)>30*2*pi) = [];
k = k*prod(abs(z(abs(z)>30*2*pi))); z(abs(z)>30*2*pi) = [];
Kc = -ZPKmod(z,p,k);
figure; bode(K); hold on; bode(Kc); 

% check whether the closed-loop specs are still satisfied, also plot the
% nonparametric loops
K.empty(); K.add(Kc); 
P.add(G{1:4});
figure(); bodemag(CL(S)); hold on; bodemag(1/MS,'--'); bodemag(info.gamma(1)/WS,'--'); title('Sensitivity S');
figure(); bodemag(CL(U)); hold on; bodemag(1/WU,'--'); title('Input sensitivity U'); 
figure(); bodemag(CL(T)); title('Closed loop T'); 

% save controller
K = std(Kc); 
save controller K; 
