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

bw = 5:10:100;
for i=1:length(bw)
    WS = Weight.LF(bw(i),2,-80); 
    MS = Weight.DC(5);
    WU = Weight.HF(bw(i)+10,2,-40);

    obj =  WS*S;
    cstr = [MS*S <= 1; 
            WU*U <= 1];

    % solve
    [~,~,info] = CL.solve(obj,cstr,K);
end

% postprocess
controllers = cellfun(@(x) shiftlow(x,9), K.content, 'un', 0); 
controllers = cellfun(@(x) removehigh(x,80*2*pi), controllers, 'un', 0); 
figure; bode(controllers{:}); 

% check whether the closed-loop specs are still satisfied
K.empty(); K.add(controllers{:}); 
%P.add(G{1:4}); % way too much info on one plot
figure(); bodemag(CL(S));  title('Sensitivity S'); hold on; legend off; %bodemag(1/MS,'--'); bodemag(info.gamma(1)/WS,'--'); 
figure(); bodemag(CL(U)); title('Input sensitivity U'); hold on; legend off; %bodemag(1/WU,'--'); 
figure(); bodemag(CL(T)); title('Closed loop T'); legend off; 

% save controller
K = cellfun(@std, controllers, 'un', 0); 
save controllers K; 
