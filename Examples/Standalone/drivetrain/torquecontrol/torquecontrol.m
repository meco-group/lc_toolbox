%% Torque controller design for mechatronic drive train
%  Laurens Jacobs, 30 Oct 2019

% initialize
clear; close all; clc; 
addpath(genpath('lc_toolbox'));

%% identification

% nonparametric data
load('frequencyRange.mat'); 
Ga = [];
G = {};
for i=[13 14 16:1:23]
    load(['rawFRF211019run' num2str(i) '.mat']);
    G = [G {FRDmod(FRF,f,'FrequencyUnit','Hz')}];
end
Gav = 1/10*(G{1}+G{2}+G{3}+G{4}+G{5}+G{6}+G{7}+G{8}+G{9}+G{10});
figure; bode(G{:}); hold on; bode(Gav,'k'); title('Measurements & average'); 

% fit
FRFW = squeeze(ones(size(Gav.ResponseData))); 
FRFW(11:500) = 0;       % don't fit to 50 Hz
FRFW(501:1000) = 1;     % fit AR/R pair
FRFW(1001:2000) = 10;   % more weight to antiresonance
FRFW(2001:end) = 0;     % don't fit from 200 Hz
settings = struct('denh',2,'denl',0,'numh',4,'numl',0,'FRFW',FRFW);
Gfit = param_ident('data',Gav,'method','nllsfdi','settings',settings);
figure; bode(Gav); hold on; bode(Gfit); title('Model'); 

% Note: this fit is not a good approximation at low frequencies (influence
% velocity controller of other motor) nor at high frequencies (due to
% noise). At low frequencies, this is not a problem as we will require
% tight control, at high frequencies we will enforce roll-off on the
% controller to ensure robustness. 

%% control design

% control configuration
P = IOSystem(Gfit);         % plant
K = IOSystem(1,1);          % controller

r = Signal(1);              % torque reference
u = P.in;                   % current command
y = P.out;                  % measured torque
e = r-y;                    % tracking error
conn = [K.in == e; 
        K.out == u];    
CL = IOSystem(P,K,conn);    % define the control configuration

% define channels
S = Channel(e/r,'Sensitivity'); 
T = Channel(y/r,'Complementary sensitivity');
U = Channel(u/r,'Input sensitivity');

% define specifications
WS = Weight.LF(100,1,-40); 
MS = Weight.DC(6);
WU = Weight.HF(100,1,-40);

obj =  WS*S;
cstr = [MS*S <= 1 ; WU*U <= 1];
    
% solve
[~,~,info] = CL.solve(obj,cstr,K);

% plot
 showall(info); % plots the channels that were in the optimization problem
 figure; bode(CL(T)); title('Closed loop'); % plot the closed-loop
 figure; bode(K); title('Controller'); % plot the controller

% postprocess
C = shiftlow(K,2*2*pi);
C = removehigh(C, 120*2*pi);
K.empty(); 
K.add(fromstd(C)); 

% check whether the closed-loop specs are still satisfied, also plot the
% nonparametric loops
P.add(G{:});
figure(); bodemag(CL(S)); hold on; bodemag(1/MS,'--'); bodemag(info.gamma(1)/WS,'--'); title('Sensitivity S');
figure(); bodemag(CL(U)); hold on; bodemag(1/WU,'--'); title('Input sensitivity U'); 
figure(); bodemag(CL(T)); title('Closed loop T'); 

% save controller
K = std(K); 
save controller K; 
